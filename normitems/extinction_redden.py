import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from numba import vectorize
@vectorize
#def A_wave(v, Rv, EBV):
def A_wave(v, Rv):
	"""
	v : the wavelength of filter, in unit "10^(-10)m"
	EBV : E(B-V)
	Rv : 3.1 for Milky Way
	"""
	x = 10**4/v
	'''
	a = np.zeros(len(x), dtype = np.float)
	b = np.zeros(len(x), dtype = np.float)
	'''
	if 0.3 <= x <= 1.1:
		"""
		infrared
		"""
		a = 0.574*x**1.61
		b = -0.527*x**1.61

	elif 1.1 <= x <= 3.3:
		"""
		Optical/NIR
		"""
		y = x - 1.82
		a = (1. + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 
			+ 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7)
		b = (1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4
			- 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7)

	elif 3.3 <= x <= 8:
		"""
		Ultraviolet
		"""
		if x<5.9:
			Fa = 0.
			Fb = 0.
		elif 5.9 <= x <= 8.:
			Fa = -0.04473*(x - 5.9)**2 - 0.009779*(x - 5.9)**3
			Fb = 0.2130*(x - 5.9)**2 + 0.1207*(x - 5.9)**3
			
		a = 1.752 - 0.316*x - 0.104/((x - 4.67)**2 + 0.341) + Fa
		b = -3.090 + 1.825*x + 1.206/((x - 4.62)**2 + 0.263) + Fb

	elif 8. <= x <= 10.:
		"""
		Far-UV
		"""
		a = -1.073 - 0.628*(x - 8) + 0.137*(x - 8)**2 - 0.070*(x - 8)**3
		b = 13.670 + 4.257*(x - 8) - 0.420*(x - 8)**2 + 0.374*(x - 8)**3

	Aw2Av = a + b/Rv
	return Aw2Av

def fig():
	w_sdss = np.linspace(3000, 11000, 100)	
	w_ref = 1/np.linspace(1, 8.6, 100)
	w_ref = w_ref*10**4

	f_sdss = A_wave(w_sdss, 3.1)
	f_ref = A_wave(w_ref, 3.1)
	x_sdss = 10**4/w_sdss
	x_ref = 10**4/w_ref

	plt.figure()
	plt.plot(x_sdss, f_sdss, ls = '--', c = 'r', label = '$SDSS$', alpha = 0.5)
	plt.plot(x_ref, f_ref, ls = '-', c = 'b', label = '$Carllide_{1989}$', alpha = 0.5)
	plt.legend(loc = 2)
	plt.xlabel('$1/ \lambda [\mu m ^{-1}]$')
	plt.ylabel('$A_{\lambda}/A_{v}$')
	plt.savefig('extinction_law.png', dpi = 600)
	plt.close()

	# compare with SFD 1998
	w_l = np.array([3546, 4925, 6335, 7799, 9294])
	band = ['u', 'g', 'r', 'i', 'z']
	com_sfd = np.array([1.579, 1.161, 0.843, 0.639, 0.453])
	com_sdss = A_wave(w_l, 3.1)
	delta = com_sdss - com_sfd

	plt.figure()
	plt.plot(w_l, com_sfd, 'ro', label = '$SFD_{1998}$', alpha = 0.5)
	plt.plot(w_l, com_sdss, 'bs', label = '$mine$', alpha = 0.5)
	plt.legend(loc = 1)
	for k in range(5):
		plt.text(w_l[k], com_sfd[k] + 0.05, s = '%.3f' % delta[k])
	plt.text(4e3, 0.45, s = '$\longrightarrow$')
	plt.text(4e3, 0.4, s = '$u \; g \; r \; i \; z$')
	plt.xlabel(r'$wavelength[\mu m]$')
	plt.ylabel(r'$A(\lambda)/A(v)$')
	plt.xlim(3.45e3, 1e4)
	plt.ylim(0.37, 1.7)
	plt.savefig('extinction_compare.png', dpi = 600)
	plt.close()

def main():
	fig()

if __name__ == "__main__":
	main()