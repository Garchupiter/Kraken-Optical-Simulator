#!/usr/bin/env Python3
from .ciddorModel import Observatory
from .refractionModels import refraction
from .dispersionModels import dispersion
from .neoslalib import sla_refro
from .observatories import observatories
from .refractivityModels import Edlen1953
import numpy as np

'''
@ Author: 		Joost van den Born
@ Contact: 		born@astron.nl
@ Description: 	Some functions for quick testing.
'''
def get_conditions(location):
	'''
		Returns the mean conditions at various site. Currently
		supported locations:
		- 'STANDARD'
		- 'CERRO_ARMAZONES'
		- 'CERRO_PARANAL'
		- 'LA_PALMA'
		- 'LA_SILLA'
		- 'LAS_CAMPANAS'
		- 'MAUNA_KEA'

		Parameters
		----------
		location 	: string
			String denoting the desired location

		Returns
		-------
		observatory : dict
			Dictionary containing the atmospheric parameters
	'''
	observatory = observatories[location]
	return observatory


def Filippenko1982(l1, l2, zenith, conditions='STANDARD', T=288.15, p=101325,
							RH=0.0, xc=450, lat=0, h=0, f=None):
	'''
		Calculates the refraction following the 1982 paper by 
		Alexei V. Filippenko.

		The refraction is based on a flat atmosphere and the refractivity
		values are calculated from a modified Edlen (1953) equation.

		Parameters
		----------
		l1 	: 	float
			Shortest wavelength of interest in micron.
		l2 	: 	float
			Longest wavelength of interest in micron.
		zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
		conditions 	: string
			Refers to the standard conditions of the 
			observatory. See get_conditions() for more info.
			If set to 'None', the conditions specified 
			after will be used.
		T 	: float (optional)
			Temperature in Kelvin
		p 	: float (optional)
			Atmospheric pressure in Pa
		RH 	: float (optional)
			Relative humidity
		xc 	: float (optional)
			CO2 density in parts per million
		lat : float (optional)
			Latitude of the observer in degrees
		h 	: float (optional)
			Altitude of the observer in meters
		f 	: float
			Water vapour pressure (default=None). If None, it is
			calculated automatically.

		Returns
		-------
		R : float
			Atmospheric refraction in degrees
	'''

	if conditions is not None:
		T, p, RH, xc, lat, h, _ = get_conditions(conditions).values()

	n1 = Edlen1953(l1, T=T, p=p, RH=RH, f=f)
	n2 = Edlen1953(l2, T=T, p=p, RH=RH, f=f)
	dR  = (n1 - n2) * np.tan(np.deg2rad(zenith))
	return np.degrees(dR)


def quick_refractive_index(l, conditions='STANDARD', T=288.15, p=101325, 
							RH=0.0, xc=450, lat=0, h=0):
	'''
		Function to quickly calculate the refractive index
		of atmospheric air at reference conditions, given 
		a wavelength l.

		Parameters
		----------
		l 	: float
			Wavelength in microns
		conditions: string
			Refers to the standard conditions of the 
			observatory. See get_conditions() for more info.
			If set to 'None', the conditions specified 
			after will be used.
		T 	: float (optional)
			Temperature in Kelvin
		p 	: float (optional)
			Atmospheric pressure in Pa
		RH 	: float (optional)
			Relative humidity
		xc 	: float (optional)
			CO2 density in parts per million
		lat : float (optional)
			Latitude of the observer in degrees
		h 	: float (optional)
			Altitude of the observer in meters
	
		Returns
		-------
		n 	: 	float
			Refractive index for the given conditions.
	'''
	if conditions is not None:
		T, p, RH, xc, lat, h, _ = get_conditions(conditions).values()

	# Initializing dispersion model
	at  = Observatory()

	# Calculating indices of refraction for l
	n 	= at.n_tph(l=l, T=T, p=p, RH=RH, xc=xc)
	return n


def quick_refraction(l, zenith, conditions='STANDARD', T=288.15, p=101325, 
							RH=0.0, xc=450, lat=0, h=0):
	'''
		Function to quickly calculate the refraction
		of atmospheric air at reference conditions, 
		given a wavelength l.

		Parameters
		----------
		l  	: 	float 
			Wavelength in microns
		zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
		conditions 	: string
			Refers to the standard conditions of the 
			observatory. See get_conditions() for more info.
			If set to 'None', the conditions specified 
			after will be used.
		T  	: 	float (optional)
			Temperature in Kelvin
		p  	: 	float (optional)
			Atmospheric pressure in Pa
		RH  	: 	float (optional)
			Relative humidity
		xc  : 	float (optional)
			CO2 density in parts per million
		lat : 	float (optional)
			Latitude of the observer in degrees
		h  	: 	float (optional)
			Altitude of the observer in meters
	
		Returns
		-------
		R 	: 	float
			Refraction in degrees
	'''
	if conditions is not None:
		T, p, RH, xc, lat, h, _ = get_conditions(conditions).values()

	# Initializing dispersion model
	at  = Observatory()

	# Calculating indices of refraction for l
	n 	= at.n_tph(l=l, T=T, p=p, RH=RH, xc=xc)

	# Density of the atmosphere (following CIPM-81/91 equations)
	rho = at.rho(p=p, T=T, RH=RH, xc=xc)

	# Initializing refraction model and setting the reduced height
	ref = refraction(lat, h)
	ref.setReducedHeight(p, rho)
	return ref.cassini(n, zenith)


def quick_dispersion(l1, l2, zenith, conditions='STANDARD', T=288.15, 
						p=101325, RH=0.0, xc=450, lat=0, h=0):
	'''
		Function to quickly calculate the dispersion
		of atmospheric air at reference conditions, 
		given a wavelength l.

		Parameters
		----------
		l1 	: 	float
			Shortest wavelength of interest in micron.
		l2 	: 	float
			Longest wavelength of interest in micron.
		zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
		conditions 	: string
			Refers to the standard conditions of the 
			observatory. See get_conditions() for more info.
			If set to 'None', the conditions specified 
			after will be used.
		T 	: 	float (optional)
			Temperature in Kelvin
		p 	: 	float (optional)
			Atmospheric pressure in Pa
		RH 	: 	float (optional)
			Relative humidity
		xc 	: 	float (optional)
			CO2 density in parts per million
		lat : 	float (optional)
			Latitude of the observer in degrees
		h 	: 	float (optional)
			Altitude of the observer in meters
	
		Returns
		-------
		dR 	: float
		 	Atmospheric dispersion in degrees
	'''
	if conditions is not None:
		T, p, RH, xc, lat, h, _ = get_conditions(conditions).values()

	# Initializing dispersion model
	at  = Observatory()

	# Calculating indices of refraction for l1 and l2
	n1 	= at.n_tph(l=l1, T=T, p=p, RH=RH, xc=xc)
	n2 	= at.n_tph(l=l2, T=T, p=p, RH=RH, xc=xc)

	# Density of the atmosphere (following CIPM-81/91 equations)
	rho = at.rho(p=p, T=T, RH=RH, xc=xc)

	# Initializing refraction model and setting the reduced height
	disp = dispersion(lat, h)
	disp.setReducedHeight(p, rho)
	return disp.cassini(n1, n2, zenith)


def slalib_refraction(l, zenith, conditions='STANDARD', T=288.15, p=101325, 
				RH=0.0, xc=450, lat=0, h=0, tlr=0.0065, eps=1e-8, use_pyslalib=False):
	'''
		Calculates the atmospheric refraction assuming the SLALIB
		package. SLALIB uses the same atmospheric model as ZEMAX
		based on the works of Hohenkerk & Sinclair and Seidelmann.
		They take the B&S refraction formula as a basis.
		
		The pySLALIB FORTRAN has been ported, so pyslalib is not
		needed, but can still be used. There are minor numerical
		differences, likely due to precision or declaration issues.

		Parameters
		---------- 
		l 	:	float
			Wavelength in microns
		zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
		conditions 	: string
			Refers to the standard conditions of the 
			observatory. See get_conditions() for more info.
			If set to 'None', the conditions specified 
			after will be used.
		T 	: 	float (optional)
			Temperature in Kelvin
		p 	: 	float (optional)
			Atmospheric pressure in Pa
		RH 	: 	float (optional)
			Relative humidity
		xc 	: 	float (optional)
			CO2 density in parts per million
		lat : 	float (optional)
			Latitude of the observer in degrees (default=0)
		h 	: 	float (optional)
			Altitude of the observer in meters
		tlr :	float
			Temperature lapse rate (default=0.0065)
		eps :	float
			Precision at which the calculation is aborted (default=1e-8)
		use_pyslalib : bool
			Flag that indicates that the user wants to use the pySlalib
			package. If False, the internal python port of the original
			FORTRAN code will be used. (default=False).
			
		Returns
		-------
		R 	: 	float
			Refraction in degrees
	'''

	if conditions is not None:
		T, p, RH, xc, lat, h, _ = get_conditions(conditions).values()

	# Converting units to the units that sla_refro() requires.
	_zenith = np.deg2rad(zenith) 	# From degrees to radians
	_P  	= p / 100 				# From Pa to mbar
	_phi  	= np.deg2rad(lat) 		# From degrees to radians

	R = []
	if type(l) == float or type(l) == int:
		l = [l]
	if use_pyslalib:
		try:
		    from pyslalib import slalib
		except ImportError as e:
		    print('pySLALIB package not found! Try installing it first.')
		for wl in l:
			# If using original slalib (e.g. PySLALIB)
			Refrac = slalib.sla_refro(_zenith, h, T, _P, RH, wl, phi, tlr, eps)
			R.append(np.degrees(Refrac))
	
	# If using internal slalib port
	else:
		for wl in l:
			Refrac = sla_refro(_zenith, h, T, _P, RH, wl, phi, tlr, eps)
			R.append(np.degrees(Refrac))

	if len(l) == 1:
		return R[0]
	else:
		return np.asarray(R)

def slalib_dispersion(l1, l2, zenith, conditions='STANDARD', T=288.15, p=101325, 
				RH=0.0, xc=450, lat=0, h=0, tlr=0.0065, eps=1e-8, use_pyslalib=False):
	'''
		Calculates the atmospheric refraction assuming the SLALIB
		package. SLALIB uses the same atmospheric model as ZEMAX
		based on the works of Hohenkerk & Sinclair and Seidelmann.
		They take the B&S refraction formula as a basis.
		
		The pySLALIB FORTRAN has been ported, so pyslalib is not
		needed, but can still be used. There are minor numerical
		differences, likely due to precision or declaration issues.

		Parameters
		---------- 
		l1 	: 	float
			Shortest wavelength of interest in micron.
		l2 	: 	float
			Longest wavelength of interest in micron.
		zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
		conditions 	: string
			Refers to the standard conditions of the 
			observatory. See get_conditions() for more info.
			If set to 'None', the conditions specified 
			after will be used.
		T 	: 	float (optional)
			Temperature in Kelvin
		p 	: 	float (optional)
			Atmospheric pressure in Pa
		RH 	: 	float (optional)
			Relative humidity
		xc 	: 	float (optional)
			CO2 density in parts per million
		lat : 	float (optional)
			Latitude of the observer in degrees (default=0)
		h 	: 	float (optional)
			Altitude of the observer in meters
		tlr :	float
			Temperature lapse rate (default=0.0065)
		eps :	float
			Precision at which the calculation is aborted (default=1e-8)
		use_pyslalib : bool
			Flag that indicates that the user wants to use the pySlalib
			package. If False, the internal python port of the original
			FORTRAN code will be used. (default=False).
			
		Returns
		-------
		dR 	: 	float
			Refraction in degrees
		'''
	R1 = slalib_refraction(l1, zenith, conditions, T, p, RH, xc,\
	 lat, h, tlr, eps, use_pyslalib)

	R2 = slalib_refraction(l2, zenith, conditions, T, p, RH, xc,\
	 lat, h, tlr, eps, use_pyslalib)

	dR = R1 - R2
	return dR
