import numpy as np 
from struct import unpack
from os import fstat
from astropy import constants as const
import math
from astropy.cosmology import Planck15 as cosmo


global pi
pi = 3.1415926

#------------------------------------------------------------------------------------

#read the number of particles
def read_npart(filename):
	file = open(filename, 'rb')	
	space_8 = unpack('<i',file.read(4))
	bname = file.read(4).decode('ascii')
	bsize = unpack('<i',file.read(4))
	space8 = file.read(8)
	npart = unpack('<iiiiii', file.read(4 * 6))
	file.close()
	return npart

#read the datas of block(except 'MASS')
def read_datas(filename,block, dtype, column, dt, file):
	npart = read_npart(filename)
	if block == 'AGE ':
		data = file.read(npart[4] * dt.itemsize * column)
	elif block == 'ZTOT' or block == 'Z   ' or block == 'ZS  ':
		if dtype == 0:
			data = file.read(npart[0] * 4 * column)
		else:
			file.seek(file.tell() + npart[0] * 4 * column)
			data = file.read(npart[4] * 4 * column)
	else:	
		if dtype == 6:
			data = file.read(np.sum(npart[:dtype]) * dt.itemsize * column)
			nsum = np.sum(npart[:])
			arr = np.ndarray(shape = (np.int32(nsum), column), dtype = dt, buffer = data)
		else:
			if dtype == 0:
				data =  file.read(npart[0] * dt.itemsize * column)
			else:
				file.seek(file.tell() + np.sum(npart[:dtype]) * dt.itemsize * column)
				data = file.read(npart[dtype] * dt.itemsize * column)

	if column != 1 and dtype != 6:
		arr = np.ndarray(shape = (np.int32(npart[dtype]), column), dtype = dt, buffer = data)
	elif column == 1 and dtype != 6:
		arr = np.ndarray(shape = (np.int32(npart[dtype])), dtype = dt, buffer = data)
	
	return arr

#find the location of block
def read_all(filename, block, dtype, column ,dt): 
	file = open(filename, 'rb')
	allsize = fstat(file.fileno()).st_size
	
	size0 = 24 + 256
	file.seek(size0)
	while file.tell() < allsize:
		space4 = file.read(4)
		bname = file.read(4).decode('ascii')
		bsize = unpack('<i',file.read(4))[0] 
		space8 = file.read(8)
		if (block != bname):
			size0 =  bsize + size0 +16
			file.seek(size0)
		else:
			return read_datas(filename,block, dtype, column, dt, file)
		pass
	file.close()
	return None

#the mass of particles
def mass_type(filename,ntype):
	npart = read_npart(filename)
	file = open(filename, 'rb')
	file.seek(20 + 24)
	massarr = unpack('<dddddd', file.read(8 * 6))

	if massarr[ntype] != 0:
		mas = massarr[ntype]
		mass = np.ones(npart[ntype])
		mass = mass * mas
	elif npart[ntype] == 0:
		mass = 0 
	else:
		size1 = 280
		for i in range(3):
			file.seek(size1)
			file.read(4)
			bname = file.read(4).decode('ascii')
			bsize = unpack('<i',file.read(4))[0] 
			size1 = size1 + bsize + 16
		locntype = 0
		for j in range(ntype):
			if massarr[j] == 0:
				locntype = locntype + npart[j] * 4		
		file.seek(size1 + 20 + locntype)
		data = file.read(npart[ntype] * 4)
		mass = np.ndarray(shape = np.int32(npart[ntype]), dtype = 'float32', buffer = data)
		
	file.close()
	return mass

#the temperture of gas 
def caltemp(filename):
	file = open(filename, 'rb')
	file.read(20 + 24 + 48)
	time = unpack('<d',file.read(8))
	temp = read_all(filename,'U   ',0,1,np.dtype('float32'))
	xH = 0.76
	yhelium = (1. - xH) / (4 * xH)
	NE = read_all(filename,'NE  ',0,1,np.dtype('float32'))
	mean_mol_weight = (1. + 4. * yhelium) / (1. + yhelium + NE)
	v_unit = 1.0e5 * np.sqrt(time)
	prtn = 1.67373522381e-24
	bk = 1.3806488e-16
	temp = temp * (5. / 3 - 1) * v_unit**2 * prtn * mean_mol_weight / bk

	temp = temp * const.k_B.value / const.e.value * 1e-3
	file.close()
	return temp

#the stars age
def calage(filename):
	a = read_all(filename,'AGE ', 4, 1, np.dtype('float32'))
	z = 1. / a - 1
	age = z
	age = cosmo.age(z).value
	return age




#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#find particles within virial radius
def distribution(xc,yc,zc,rv,filename,nn):# h size of interval

	#read datas
	npart = read_npart(filename)
	pos0 = read_all(filename,'POS ', 0, 3, np.dtype('float32'))
	#pos4 = read_all(filename,'POS ', 4, 3, np.dtype('float32'))

	# sfr = read_all(filename,'SFR ',0,1,np.dtype('float32'))
	temp = caltemp(filename)

	# age = calage(filename)
	# age = cosmo.age(0).value - read_all(filename,'AGE ', 4, 1, np.dtype('float32')) / 1e9 # the age for Gadget3PESPH
	# np.savetxt('/home/qyli/yearG.txt',age)
	

	# metalstars = read_all(filename,'Z   ',4,1,np.dtype('float32'))
	metalgas = read_all(filename,'Z   ',0,1,np.dtype('float32'))
	# metalstars = np.sum(read_all(filename,'ZS  ',4,4,np.dtype('float32')), axis = 1) # the metallicity for Gadget3PESPH
	# metalgas = np.sum(read_all(filename,'ZS  ',0,4,np.dtype('float32')), axis = 1) # the metallicity for Gadget3PESPH

	mass_0 = mass_type(filename,0)
	# mass_4 = mass_type(filename,4)
	

	h = rv / nn

	mass0 = np.zeros(nn)
	temp_sum = np.zeros(nn)
	age_sum = np.zeros(nn)
	metalgas_sum = np.zeros(nn)

	rr  = np.zeros(nn)
	count0 = np.zeros(nn)

	i = 0
	while npart[0] > i:
		rc = math.sqrt( (pos0[i][0] - xc)**2 + (pos0[i][1] - yc)**2 + (pos0[i][2] - zc)**2) 
		
		if rc < rv:
			n = int(rc / h)
			mass0[n] = mass_0[i] + mass0[n]
			temp_sum[n] = temp[i]*mass_0[i] + temp_sum[n] # mass-weighted temperature  
			metalgas_sum[n] = metalgas[i]*mass_0[i] + metalgas_sum[n] # mass-weighted metallicity

			rr[n] = rc / rv * mass_0[i] + rr[n]
			count0[n] += 1
		i += 1

	temp00 = temp_sum / mass0
	metalgas00 = metalgas_sum / mass0
	rr = rr / mass0

	for k in range(nn):
		dv = 4. / 3 * math.pi *(pow((k+1)*h,3) - pow(k*h,3))
		mass0[k] = mass0[k] / dv

	return rr,count0,mass0,temp00,metalgas00

	
#-------------------------------------------------------------------------------------

def main(nn):
	
	speciname = '/home/nifty2014/TheThreeHundred/catalogues/AHF/complete_sample/G3X_Mass_snap_128-center-cluster.txt'
	data = np.loadtxt(speciname)
	size = data.shape[0]

	density0 = np.zeros(shape = (324,nn))
	temp = np.zeros(shape = (324,nn))
	gasz = np.zeros(shape = (324,nn))
	rr0 = np.zeros(shape = (324,nn))
	ngas = np.zeros(shape = (324,nn))

	ii = 0
	while  ii < size:
		region = int(data[ii][0])

		xc = data[ii][3]
		yc = data[ii][4]
		zc = data[ii][5]
		rv = data[ii][10]


		if region <= 9:
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_000%s/snap_128' %region
		if  9 < region <= 99:
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_00%s/snap_128' %region
		if   region > 99:
			filename = '/home/nifty2014/TheThreeHundred/simulation/GadgetX/NewMDCLUSTER_0%s/snap_128' %region

	
		rr,count0,mass0,temp00,metalgas00 = distribution(xc,yc,zc,rv,filename,nn)

		density0[ii] = mass0
		temp[ii] = temp00
		gasz[ii] = metalgas00
		rr0[ii] = rr
		ngas[ii] = count0

		ii += 1
		print(ii)

	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/GX500_snap_density0.txt', density0)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/GX500_snap_temp.txt', temp)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/GX500_snap_gasz.txt', gasz)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/GX500_snap_rr0.txt', rr0)
	np.savetxt('/home/nifty2014/TheThreeHundred/playground/qingyang/simdata/GX500_snap_ngas.txt', ngas)
		
		
main(20)







