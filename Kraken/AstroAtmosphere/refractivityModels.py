#!/usr/bin/env Python3
import numpy as np

'''
@ Author:		Joost van den Born 
@ Contact:		born@astron.nl
@ Description:	A comparison of various atmospheric refraction models.
				Also these are used to calculate the difference in 
				atmospheric dispersion, given a set of conditions.
'''

def BarrellAndSears(l, T=288.15, p=101325, RH=0, f=None):
	'''
		From Barrell & Sears, 1939
		Equation 7.7, page 52
		For atmospheric air (i.e. CO2 = 300 ppm)
		This equation is applicable for the  ranges T=10-30C and p=720-800 mmHg.

		Parameters
		---------- 
		l 	:	float
			Wavelength in microns
		T 	: 	float
			Temperature in Kelvin
		p 	: 	float
			Atmospheric pressure in Pascal
		RH 	: 	float
			Relative humidity (0=<H=<1)
		f 	: 	float
			Water vapour pressure in Pascal

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	if f is None:
		f = RH * svp(T)
	f_mmHg = hPa2mmHg(f/100)
	T = K2C(T)
	P = hPa2mmHg(p/100)
	R = (0.378125 + 0.0021414 / l**2 + 0.00001793 / l**4) * (P*(1 + (1.049 - 0.0157*T) * P * 1e-6)/(1 + 0.003661*T)) - (0.0624 - 0.000680 / l**2) * f_mmHg / (1 + 0.003661*T)
	n = R/1e6 + 1
	return n


def Ciddor(l, CO2=450):
	'''
		From Ciddor, 1996
		Equations 1 and 2
		At default conditions only. For more flexibility use
		Observatory.n_tph()

		Parameters
		---------- 
		l 	: 	float
			Wavelength in microns
		CO2 : 	float
			CO2 density in parts per million

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	s = 1 / l
	R = 5792105 / (238.0185 - s**2) + 167917 / (57.362 - s**2)
	if CO2 != 450:
		R = R * (1 + 0.534e-6 * (CO2 - 450))  
	n = R/1e8 + 1
	return n


def Edlen1953(l, T=288.15, p=101325, RH=0, f=None):
	'''
		From Edlen, 1953
		Equation 2

		At standard conditions only.
		Following Filippenko (1982), eqs. 2 & 3 we also include
		the	modification proposed by Barrell in 1951.

		Parameters
		---------- 
		l 	: 	float
			Wavelength in microns
		T 	: 	float
			Temperature in Kelvin
		p 	: 	float
			Atmospheric pressure in Pascal
		RH 	: 	float
			Relative humidity (0=<H=<1)
		f 	: 	float
			Water vapor pressure in Pascal

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	s = 1 / l
	R = 6432.8 + 2949810 / (146 - s**2) + 25540/(41-s**2)
	n_standard = R/1e8 + 1

	pHg = hPa2mmHg(p/100)
	Tc  = K2C(T)
	if f is None:
		f = RH * svp(T)
	f_mmHg = hPa2mmHg(f/100)

	# Changes in temperature and pressure, from Barrell (1951)
	n_tp = R/1e8 * (pHg * (1 + (1.049 - 0.0157*Tc) * 1e-6 * pHg ) ) / ( 720.883 * (1 + 0.003661 * Tc)) + 1

	# Correction for water vapour pressure
	n_svp = (n_tp - 1) * 1e6 - ( 0.0624 - 0.000680 * s**2 ) / (1 + 0.003661 * Tc) * f_mmHg
	n = n_svp/1e6 + 1
	return n


def Edlen1966(l, T=288.15, p=101325, RH=0, CO2=300, f=None):
	'''
		From Edlen, 1966
		Equations 1, 12, 17, 22

		At standard conditions only

		Parameters
		---------- 
		l 	: 	float
			Wavelength in microns
		T 	: 	float
			Temperature in Kelvin
		p 	: 	float
			Atmospheric pressure in Pascal
		RH 	: 	float
			Relative humidity (0=<H=<1)
		CO2 : 	float
			CO2 density in parts per million
		f 	: 	float
			Water vapor pressure in Pascal

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	s = 1 / l
	R = 8342.13 + 2406030 / (130 - s**2) + 15997 / (38.9 - s**2)
	n = R/1e8 + 1

	# Influence of temperature and pressure, eq. 12
	Tc = K2C(T)
	Ptorr = Pa2Torr(p)
	R_tp = ( Ptorr * R/1e8 ) / 720.775 * (1 + Ptorr * (0.817 - 0.0133 * Tc) * 1e-6) / ( 1 + 0.0036610 * Tc)

	# Allow for a CO2 fraction, eq. 17, but we include T and P dependence.
	# This is okay from "Asthe refractivity of CO, is about 50% higher than
	# that of air this CO,-content increases the refractivity over that of 
	# CO2-free air by approximately 0.5 x 0.0003 x (n-1)s" and "As the 
	# density factor of CO, is not much different from that of air this ratio r
	# can be taken as independent of t and p ."
	R_x = R_tp * (1 + 0.540 * (CO2/1e6 - 0.0003))
	n_x = R_x + 1

	# Correction for water vapour, eq. 22
	if f is None:
		f = RH * svp(T)
	ftorr = Pa2Torr(f)

	n = n_x - ftorr * (5.7224 - 0.0457 * s**2) * 1e-8
	return n


def Owens(l, T=288.15, p=101325, RH=0):
	'''
		From Owens, 1967
		Equation 32
		Note that the pressure P and H are in torricelli (Torr),
		however, this is almost precisely 1 Torr = 1 mmHg

		Parameters
		----------
		l 	: 	float
			wavelength in micron
		T 	: 	float
			Temperature in Kelvin
		p 	: 	float
			Atmospheric pressure in Pa
		RH 	: 	float
			Relative humdity (0<=H<=1)

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	P = hPa2mmHg(p/100)
	T = K2C(T)
	s = 1 / l
	R = 8342.13 + 2406030 / (130 - s**2) + 15997 / (38.9 - s**2) * (P/720.775) * ((1 + P * (0.817 - 0.0133) * 1e-6)/(1 + 0.0036610*T)) - RH * (5.722 - 0.0457*s**2)
	n = R/1e8 + 1
	return n


def BonschPotulski(l, T=288.15, p=101325, RH=0, CO2=400, f=None):
	'''
		From Bonsch & Potulski, 1998
		Equations 6, 7, 8 and 9a
		
		Parameters
		---------- 
		l 	: 	float
			Wavelength in microns
		T 	: 	float
			Temperature in Kelvin
		p 	: 	float
			Atmospheric pressure in Pa
		RH 	: 	float
			Relative humidity (0<=H<=1)
		CO2 : 	float
			Particles per million of CO2
		f 	: 	float
			Water vapor pressure in Pascal

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	x = CO2 / 1e6
	Tc = K2C(T) 
	s = 1 / l

	Rs = 8091.37 + 2333983 / (130 - s**2) + 15518 / (38.9 - s**2)
	Rs = Rs / 1e8
	Rx = Rs * (1 + 0.5327 * (x - 0.0004))
	Rtp = Rx * p / 93214.60 * (1 + 1e-8 * (0.5953 - 0.009876 * Tc)*p) / (1 + 0.0036610 * Tc)

	if f is None:
		f = RH * svp(T)
	
	n_tp = Rtp + 1
	n = n_tp - f * (3.8020 - 0.0384 * s**2) * 1e-10
	return n


def BirchDowns(l):
	'''
		From Birch & Downs, 1993 and 1994 (correction)
		Equations 1 and 2

		Only at standard conditions

		Parameters
		----------
		l 	: 	float
			wavelength in micron

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	s = 1 / l
	R = 8342.54 + 2406147 / (130 - s**2) + 15998 / (38.9-s**2)
	n = R/1e8 + 1
	return n


def PeckReeder(l):
	'''
		From Peck & Reeder, 1972
		Equation 3

		Only at standard conditions

		Parameters
		----------
		l 	: 	float
			wavelength in micron

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	s = 1 / l
	R = 8060.51 + 2480990 / (132.274 - s**2) + 17455.7 / (39.32957 - s**2)
	n = R/1e8 + 1
	return n


def HohenkerkAndSinclair(l, T=288.15, p=101325, RH=0, h=0, lat=0):
	'''
		Calcalation of the refractive index following Hohenkerk & Sinclair
		(1985), which is the same one as used in SLALIB. 

		This equation is based on the Barrell & Sears equation, but with
		slightly different values. The refractive index values here seem 
		to be a bit higher than the original Barrel and Sears values.

		Parameters
		----------
		l 	:	float
			Wavelength in microns
		T 	:	float
			Temperature in Kelvin
		p 	:	float
			Atmospheric pressure in Pa
		RH 	:	float
			fractional humidity (0=<H=<1)
		h 	:	float
			Altitude of the observer (default=0 m)
		lat :	float
			Latitude of the observer (default=-24)
		
		Returns
		-------
		n 	:	float
			Refractive index of atmospheric air at the given conditions
	'''

	# Some constants
	alpha = 0.0065
	delta = 18.36
	R = 8314.36
	Md = 28.966
	Mw = 18.016

	# Pressure in millibars
	P0 = p/100 

	# Local acceleration due to gravity
	g_local = 9.784 * (1 - 0.0026 * np.cos(2*np.deg2rad(lat))) - 0.00000028*h
	gamma = g_local * Md / (R * alpha)

	# Partial pressure of water vapour at observer
	Pw = RH * (T/247.1)**18.36
	
	# Refractivity * 1e6 at reference conditions
	A = (287.604 + 1.6288*l**-2 + 0.0136*l**-4) * 273.15/1013.25
	# A = (287.5522 + 1.6288*l**-2 + 0.0136*l**-4) * 273.15/1013.25
	
	P = (P0 + (1 - Mw/Md) * gamma/(delta - gamma) * Pw )\
		- (1 - Mw/Md) * gamma / (delta-gamma) * Pw
	
	# Refractive index at given conditions
	n = 1 + 1e-6 * (A * P - 11.2684*Pw)/T
	return n


def slalib(l, T=288.15, p=101325, RH=0, h=0, lat=-24, tlr=0.0065, eps=1e-10):
	'''
		NOTE: This function still works for now, but is considered DEPRECATED.
			  Use BarrellAndSears() or HohenkerkAndSinclair() instead.

		SLALIB uses the same atmospheric model as ZEMAX
		based on the works of Hohenkerk & Sinclair and Seidelmann.
		They take the B&S refraction formula as a basis.
		Requires the pySLALIB package to work.

		Parameters
		---------- 
		l 	:	float
			Wavelength in microns
		T 	:	float
			Temperature in Kelvin
		p 	:	float
			Atmospheric pressure in Pa
		RH 	:	float
			Relative humidity (0=<H=<1)
		h 	:	float
			Altitude of the observer (default=0 m)
		lat :	float
			Latitude of the observer (default=-24)
		tlr :	float
			temperature lapse rate (default=0.0065)
		eps :	float
			precision at which the calculation is aborted (default=1e-10)

		Returns
		-------
		n 	:	float
			Refractive index of atmospheric air at the given conditions
	
		Notes
		-----
		1. 	I don't understand why this function works. It calculates the 
			refraction and somehow returns the refractive index. Perhaps 
			this is a very particular set of input parameters that make it 
			work. It does closely follow BarrellAndSears() for standard
			conditions.
	'''
	try:
	    from pyslalib import slalib
	except ImportError as e:
	    print('pySLALIB package not found! Try installing it first.\nReturning n=1.')
	    return 1
	
	zObs = np.deg2rad(45.06563)
	p = p/100
	R = []
	if type(l) == float:
		l = [l]
	for wl in l:
		R.append(slalib.sla_refro(zObs, h, T, p, RH, wl, lat, tlr, eps))
	if len(l) == 1:
		return R[0]+1
	else:
		return np.asarray(R)+1


def Mathar(l, T=288.15, p=101325, RH=0):
	'''
		Refractivity model by Mathar (2007), valid above 1.4 microns
		Code taken from
		https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/scripts/Mathar%202007%20-%20Air%201.3-2.5um.py
		(better known from refractiveindex.info). The code was subsequently 
		minimally modified to make it suitable for this package.


		Parameters
		---------- 
		l 	: 	float
			Wavelength in microns
		T 	: 	float
			Temperature in Kelvin
		p 	: 	float
			Atmospheric pressure in Pa
		H 	: 	float
			fractional humidity (0=<H=<1)

		Returns
		-------
		n 	: 	float
			Refractive index of atmospheric air at the given conditions
	'''
	if l < 1.3 or l > 24:
		print('Using the Mathar model at this wavelength (l={:.3f}) is not recommended!'.format(l))

	if l <= 2.65:
		########################################################################
		# model parameters for wavelengths between 1.3 and 2.5 microns
		cref = [ 0.200192e-3,  0.113474e-9,  -0.424595e-14,  0.100957e-16, -0.293315e-20,  0.307228e-24] # cm^j
		cT   = [ 0.588625e-1, -0.385766e-7,   0.888019e-10, -0.567650e-13,  0.166615e-16, -0.174845e-20] # cm^j · K
		cTT  = [-3.01513,      0.406167e-3,  -0.514544e-6,   0.343161e-9,  -0.101189e-12,  0.106749e-16] # cm^j · K^2
		cH   = [-0.103945e-7,  0.136858e-11, -0.171039e-14,  0.112908e-17, -0.329925e-21,  0.344747e-25] # cm^j · %^-1
		cHH  = [ 0.573256e-12, 0.186367e-16, -0.228150e-19,  0.150947e-22, -0.441214e-26,  0.461209e-30] # cm^j · %^-2
		cp   = [ 0.267085e-8,  0.135941e-14,  0.135295e-18,  0.818218e-23, -0.222957e-26,  0.249964e-30] # cm^j · Pa^-1
		cpp  = [ 0.609186e-17, 0.519024e-23, -0.419477e-27,  0.434120e-30, -0.122445e-33,  0.134816e-37] # cm^j · Pa^-2
		cTH  = [ 0.497859e-4, -0.661752e-8,   0.832034e-11, -0.551793e-14,  0.161899e-17, -0.169901e-21] # cm^j · K · %^-1
		cTp  = [ 0.779176e-6,  0.396499e-12,  0.395114e-16,  0.233587e-20, -0.636441e-24,  0.716868e-28] # cm^j · K · Pa^-1
		cHp  = [-0.206567e-15, 0.106141e-20, -0.149982e-23,  0.984046e-27, -0.288266e-30,  0.299105e-34] # cm^j · %^-1 · Pa^-1

		sref = 1e4/2.25    # cm^−1, reference wavenumber
		Tref = 273.15+17.5 # K, reference temperature
		pref = 75000       # Pa, reference pressure
		Href = 10          #%, reference humdity in percent
		########################################################################
	
	if l > 2.65 and l <= 4.275:
		########################################################################
		# model parameters for wavelengths between 2.8 and 4.2 microns
		cref = [ 0.200049e-3,  0.145221e-9,   0.250951e-12, -0.745834e-15, -0.161432e-17,  0.352780e-20] # cm^j
		cT   = [ 0.588432e-1, -0.825182e-7,   0.137982e-9,   0.352420e-13, -0.730651e-15, -0.167911e-18] # cm^j · K
		cTT  = [-3.13579,      0.694124e-3,  -0.500604e-6,  -0.116668e-8,   0.209644e-11,  0.591037e-14] # cm^j · K^2
		cH   = [-0.108142e-7,  0.230102e-11, -0.154652e-14, -0.323014e-17,  0.630616e-20,  0.173880e-22] # cm^j · %^-1
		cHH  = [ 0.586812e-12, 0.312198e-16, -0.197792e-19, -0.461945e-22,  0.788398e-25,  0.245580e-27] # cm^j · %^-2
		cp   = [ 0.266900e-8,  0.168162e-14,  0.353075e-17, -0.963455e-20, -0.223079e-22,  0.453166e-25] # cm^j · Pa^-1
		cpp  = [ 0.608860e-17, 0.461560e-22,  0.184282e-24, -0.524471e-27, -0.121299e-29,  0.246512e-32] # cm^j · Pa^-2
		cTH  = [ 0.517962e-4, -0.112149e-7,   0.776507e-11,  0.172569e-13, -0.320582e-16, -0.899435e-19] # cm^j · K · %^-1
		cTp  = [ 0.778638e-6,  0.446396e-12,  0.784600e-15, -0.195151e-17, -0.542083e-20,  0.103530e-22] # cm^j · K · Pa^-1
		cHp  = [-0.217243e-15, 0.104747e-20, -0.523689e-23,  0.817386e-26,  0.309913e-28, -0.363491e-31] # cm^j · %^-1 · Pa^-1

		sref = 1e4/3.4     # cm^−1, reference wavenumber
		Tref = 273.15+17.5 # K, reference temperature
		pref = 75000       # Pa, reference pressure
		Href = 10          #%, reference humdity in percent
		########################################################################

	if l > 4.275 and l <= 6.35:
		########################################################################
		# model parameters for wavelengths between 4.35 and 5.2 microns
		cref = [ 0.200020e-3,  0.275346e-9,   0.325702e-12, -0.693603e-14,  0.285610e-17,  0.338758e-18] # cm^j
		cT   = [ 0.590035e-1, -0.375764e-6,   0.134585e-9,   0.124316e-11,  0.508510e-13, -0.189245e-15] # cm^j · K
		cTT  = [-4.09830,      0.250037e-2,   0.275187e-6,  -0.653398e-8,  -0.310589e-9,   0.127747e-11] # cm^j · K^2
		cH   = [-0.140463e-7,  0.839350e-11, -0.190929e-14, -0.121399e-16, -0.898863e-18,  0.364662e-20] # cm^j · %^-1
		cHH  = [ 0.543605e-12, 0.112802e-15, -0.229979e-19, -0.191450e-21, -0.120352e-22,  0.500955e-25] # cm^j · %^-2
		cp   = [ 0.266898e-8,  0.273629e-14,  0.463466e-17, -0.916894e-23,  0.136685e-21,  0.413687e-23] # cm^j · Pa^-1
		cpp  = [ 0.610706e-17, 0.116620e-21,  0.244736e-24, -0.497682e-26,  0.742024e-29,  0.224625e-30] # cm^j · Pa^-2
		cTH  = [ 0.674488e-4, -0.406775e-7,   0.289063e-11,  0.819898e-13,  0.468386e-14, -0.191182e-16] # cm^j · K · %^-1
		cTp  = [ 0.778627e-6,  0.593296e-12,  0.145042e-14,  0.489815e-17,  0.327941e-19,  0.128020e-21] # cm^j · K · Pa^-1
		cHp  = [-0.211676e-15, 0.487921e-20, -0.682545e-23,  0.942802e-25, -0.946422e-27, -0.153682e-29] # cm^j · %^-1 · Pa^-1
		σref = 1e4/4.8     # cm^−1, reference wavenumber
		Tref = 273.15+17.5 # K, reference temperature
		pref = 75000       # Pa, reference pressure
		Href = 10          #%, reference humdity in percent
		########################################################################

	if l > 6.35 and l <= 15.05:
		########################################################################
		# model parameters for wavelengths between 7.5 and 14.1 microns
		cref = [ 0.199885e-3,  0.344739e-9,  -0.273714e-12,  0.393383e-15, -0.569488e-17,  0.164556e-19] # cm^j
		cT   = [ 0.593900e-1, -0.172226e-5,   0.237654e-8,  -0.381812e-11,  0.305050e-14, -0.157464e-16] # cm^j · K
		cTT  = [-6.50355,      0.103830e-1,  -0.139464e-4,   0.220077e-7,  -0.272412e-10,  0.126364e-12] # cm^j · K^2
		cH   = [-0.221938e-7,  0.347377e-10, -0.465991e-13,  0.735848e-16, -0.897119e-19,  0.380817e-21] # cm^j · %^-1
		cHH  = [ 0.393524e-12, 0.464083e-15, -0.621764e-18,  0.981126e-21, -0.121384e-23,  0.515111e-26] # cm^j · %^-2
		cp   = [ 0.266809e-8,  0.695247e-15,  0.159070e-17, -0.303451e-20, -0.661489e-22,  0.178226e-24] # cm^j · Pa^-1
		cpp  = [ 0.610508e-17, 0.227694e-22,  0.786323e-25, -0.174448e-27, -0.359791e-29,  0.978307e-32] # cm^j · Pa^-2
		cTH  = [ 0.106776e-3, -0.168516e-6,   0.226201e-9,  -0.356457e-12,  0.437980e-15, -0.194545e-17] # cm^j · K · %^-1
		cTp  = [ 0.77368e-6,   0.216404e-12,  0.581805e-15, -0.189618e-17, -0.198869e-19,  0.589381e-22] # cm^j · K · Pa^-1
		cHp  = [-0.206365e-15, 0.300234e-19, -0.426519e-22,  0.684306e-25, -0.467320e-29,  0.126117e-30] # cm^j · %^-1 · Pa^-1
		σref = 1e4/10.1    # cm^−1, reference wavenumber
		Tref = 273.15+17.5 # K, reference temperature
		pref = 75000       # Pa, reference pressure
		Href = 10          #%, reference humdity in percent
		########################################################################


	if l > 15.05:
		########################################################################
		# model parameters for wavelengths between 16 and 24 microns
		cref = [ 0.199436e-3,  0.299123e-8,  -0.214862e-10,  0.143338e-12,  0.122398e-14, -0.114628e-16] # cm^j
		# something seems to be wrong with cT...
		cT   = [ 0.621723e-1, -0.177074e-4,   0.152213e-6,  -0.954584-9,   -0.996706e-11,  0.921476e-13] # cm^j · K
		cTT  = [-23.2409,      0.108557,     -0.102439e-2,   0.634072e-5,   0.762517e-7,  -0.675587e-9 ] # cm^j · K^2
		cH   = [-0.772707e-7,  0.347237e-9,  -0.272675e-11,  0.170858e-13,  0.156889e-15, -0.150004e-17] # cm^j · %^-1
		cHH  = [-0.326604e-12, 0.463606e-14, -0.364272e-16,  0.228756e-18,  0.209502e-20, -0.200547e-22] # cm^j · %^-2
		cp   = [ 0.266827e-8,  0.120788e-14,  0.522646e-17,  0.783027e-19,  0.753235e-21, -0.228819e-24] # cm^j · Pa^-1
		cpp  = [ 0.613675e-17, 0.585494e-22,  0.286055e-24,  0.425193e-26,  0.413455e-28, -0.812941e-32] # cm^j · Pa^-2
		cTH  = [ 0.375974e-3, -0.171849e-5,   0.146704e-7,  -0.917231e-10, -0.955922e-12,  0.880502e-14] # cm^j · K · %^-1
		cTp  = [ 0.778436e-6,  0.461840e-12,  0.306229e-14, -0.623183e-16, -0.161119e-18,  0.800756e-20] # cm^j · K · Pa^-1
		cHp  = [-0.272614e-15, 0.304662e-18, -0.239590e-20,  0.149285e-22,  0.136086e-24, -0.130999e-26] # cm^j · %^-1 · Pa^-1
		σref = 1e4/20      # cm^−1, reference wavenumber
		Tref = 273.15+17.5 # K, reference temperature
		pref = 75000       # Pa, reference pressure
		Href = 10          # %, reference humdity in percent
		########################################################################

	s = 1e4/l # cm^-1
	n = 1
	H = RH*100

	for j in range(0, 6):
		n += ( cref[j] + cT[j]*(1/T-1/Tref) + cTT[j]*(1/T-1/Tref)**2
			+ cH[j]*(H-Href) + cHH[j]*(H-Href)**2
			+ cp[j]*(p-pref) + cpp[j]*(p-pref)**2
			+ cTH[j]*(1/T-1/Tref)*(H-Href)
			+ cTp[j]*(1/T-1/Tref)*(p-pref)
			+ cHp[j]*(H-Href)*(p-pref) ) * (s-sref)**j   
	return n 


# Other necessary conversion or calculations
def K2C(K):
	'''
		Kelvin to Celsius conversion

		Parameters
		----------
		K : float
			Temperature in Kelvin

		Returns
		-------
		C : float
			Temperature in Celsius
	'''
	return K - 273.15


def hPa2mmHg(hPa):
	'''
		hPa to mmHg conversion

		Parameters
		----------
		hPa 	:	float
			Atmospheric pressure in hPa (hecto Pascal)

		Returns
		-------
		mmHg 	: 	float
			Atmospheric pressure in mmHg
	'''
	return 0.7500616827 * hPa


def Pa2Torr(Pa):
	'''
		Pa to Torr conversion

		Parameters
		----------
		Pa 	:	float
			Atmospheric pressure in Pa (Pascal)

		Returns
		-------
		Torr 	: 	float
			Atmospheric pressure in Torr
	'''

	return Pa * 760/101325


def svp(T):
	'''
		Saturation vapour pressure in Pa.

		Parameters
		----------
		T 	: float
			Atmospheric temperature in Kelvin

		Returns
		-------
		svp	: float
			Saturation vapour pressure in Pa
	'''
	_A =  1.2378847e-5 # K^-2
	_B = -1.9121316e-2 # K^-2
	_C = 33.93711047
	_D = -6.3431645e3  # K^-2
	return np.exp(_A*T**2 + _B*T + _C + _D/T)