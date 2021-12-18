import h5py
import numpy as np
import pandas as pds

def h5list( f, tab ):

	default_doc = open('/home/xkchen/h5_file_header_record.txt', 'w')

	print(tab, 'Group:', f.name, 'len:%d'%len(f), file = default_doc )

	mysp2 = tab[:-1]+ '  |-*'

	for vv in f.attrs.keys():
		print( mysp2, end = ' ', file = default_doc )
		print( '%s = %s'% (vv,f.attrs[vv]), file = default_doc )

	mysp=tab[:-1] + '  |-'
	for k in f.keys():
		d = f[k]
		if isinstance(d, h5py.Group):
			h5list(d, mysp)

		elif isinstance(d,h5py.Dataset):
			print(mysp,'Dataset:',d.name,'(size:%d)'%d.size, file = default_doc )
			mysp1=mysp[:-1]+ '  |-'

			print(mysp1,'(dtype=%s)'%d.dtype, file = default_doc )

			if d.dtype.names is not None:
				print(mysp,end = ' ', file = default_doc )
				for vv in d.dtype.names:
					print(vv,end = ',', file = default_doc )
				# print() 

			mysp2=mysp1[:-1]+ '  |-*'

			for vv in d.attrs.keys():
				print(mysp2,end = ' ', file = default_doc )
				try:
					print('%s = %s'% (vv,d.attrs[vv]), file = default_doc )
				except TypeError as e:
					print('%s = %s'% (vv,e), file = default_doc )
				except:
					print('%s = ?? Other ERR'% (vv,), file = default_doc )
		else:
			print('??->',d,'Unkown Object!', file = default_doc )

	default_doc.close()

	return

home = '/media/xkchen/My Passport/data/CSST_mock/'
out_path = '/home/xkchen/mywork/CSST/test_cluster/'

## ... list header information
f = h5py.File(home + 'cluster_lensing_cat.hdf5','r')
h5list(f,'')
f.close()

