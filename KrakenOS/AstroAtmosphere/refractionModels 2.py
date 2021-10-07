#!/usr/bin/env Python3
import numpy as np
from scipy import integrate
import os
from .ciddorModel import Observatory

'''
@ Author:		Joost van den Born 
@ Contact:		born@astron.nl
@ Description:	This script offers a complete model of the atmospheric refraction
				and dispersion. It uses the Ciddor1996 paper to calculate the 
				index of refraction for a wide range of conditions, with corresponding error propagation if needed.
				The actual refraction is calculated by any of the following models:
				- Plane parallel atmosphere
				- Cassini spherical homogeneous atmosphere
				- Mathar's Barometric exponential model
				- Other methods described in Corbard et al., 2019.
'''

class refraction(Observatory):
	'''
		Class that sets the observer location.
		These location attributes are then used to calculate the 
		appropriate atmospheric refraction.

		See also:
		Observatory()
		dispersion()
		...

		Attributes
		-----------
		lat : 	float
			Latitude of the observer in degrees
		h 	: 	float
			Altitude above mean sea level of the observer in meters
		Other Observatory() attributes

		Methods
		-------
		cassini(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric refraction assuming Cassini's homegeneous
			atmosphere model.
		cassiniError(z, l1, l2, T, p, RH, xc, dl1=0, dl2=0, dT=0.2, dP=20, 
					dRH=0.02, dCO2=20, dz=0, lat=None, h=None)
			Returns the uncertainty in the atmospheric refraction, assuming 
			Cassini's homegeneous atmosphere model and the given parameters.
		corbard(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric refraction following the method described 
			in Corbard et al., 2019.
		gammaCoefficients(n1, reduced_height=None)
			Helper function for matharExponential(), that calculates the gamma
			coefficients.
		H_isotherm(T)
			Sets the reduced height for the conditions at the observer, assuming
			an isothermal atmosphere.
		matharExponential(n1, zenith, reduced_height=None)
			Returns the atmospheric refraction assuming	Mathar's barometric 
			exponential model.
		oriani(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric refraction following the classic version
			of Oriani's theorem.	
		planeParallel(n1, n2, zenith)
			Returns the atmospheric refraction assuming a plane parallel
			atmosphere.
		psi(x)
			Submethod for corbard().
		refractionIntegral(n1, n2, zenith, R0=None, useStandardAtmosphere=True, heightData=None, rhoData=None)
			Returns the atmospheric refraction following the refraction integral.
		setReducedHeight(p, rho)
			Sets the reduced atmospheric height.
		set_H(p, rho)
			Sets the reduced height for the conditions at the observer, assuming
			an adiabatic atmosphere.
		set_rc()
			Sets the local radius of curvature of the Earth at the observer.			
		tan5(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric refraction following Oriani's theorem,
			including a fifth order term in the expansion.				
		Other Observatory() methods
	'''
	def __init__(self, lat, h):
		Observatory.__init__(self)
		self.h  	= h 					# height above sea level in meters
		self.lat 	= lat					# Latitude of the observer in degrees
		self.rc 	= self.set_rc() * 1000	# Radius of curvature of the earth at the observer
		# self.H2 	= self.H_isotherm(T) 	# Reduced height of isotherm atmosphere


	def cassini(self, n1, zenith, reduced_height=None):
		'''
			Refraction of spherical atmosphere, derived from geometry using Sine law

			Parameters
			----------
			n1 : float
				Refractive index at the wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			R : float
				Atmospheric refraction in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H
		
		_R1 = np.arcsin(n1 * self.rc * np.sin(np.deg2rad(zenith)) / (self.rc + H)) - np.arcsin(self.rc * np.sin(np.deg2rad(zenith)) / (self.rc + H) )
		return np.degrees(_R1)


	def cassiniError(self, n, dn, zenith, dz=0, reduced_height=None):
		'''
		Calculates the uncertainty in the refraction angle in degrees.
		See Observatory().n_tph() and Observatory().dn_tph() on how to
		calculate the refractive index and it's uncertainty.

		Parameters
		----------
		n 	: 	float
			Refractive index of air
		dn 	: 	float
			Uncertainty in the refractive index of air
		zenith : float
			Angle of observation in degrees.
			zenith = 0 means the observation is directly overhead.
		dz  : 	float (optional)
			Uncertainty in the zenith angle, in degrees.
		reduced_height : float (optional)
			Reduced height of the atmosphere. Default value is calculated
			from the object attributes.

		Returns
		-------
		dR 	: float
			Uncertainty in the atmospheric refraction in degrees.
		'''
		if reduced_height:
			_H = reduced_height
		else:
			_H = self.H
		_zenith = np.deg2rad(zenith)
		_dz 	= np.deg2rad(dz)
		_R 		= self.rc
		_dRdn 	= (_R * np.sin(_zenith)) / ((_R + _H) * np.sqrt(1 - ((n**2 * _R**2 * np.sin(_zenith)**2)/(_R + _H)**2)))

		_dRdz 	= (_R * np.cos(_zenith)) / (_R + _H) * ( n / np.sqrt(1-(n**2 * _R**2 * np.sin(_zenith)**2)/(_R + _H)**2) - 1 / np.sqrt(1-(_R**2 * np.sin(_zenith)**2)/(_R + _H)**2) )

		_dR 	= np.sqrt( (_dRdn * dn)**2 + (_dRdz * _dz)**2)

		return np.degrees(_dR)


	def corbard(self, n1, zenith, reduced_height=None):
		'''
			Corbard et al., 2019, mention an additional formula based on the error function.
			Oriani's theorem can be derived from this equation by 'keeping only the three first
			terms of its asymptotic expansion'.

			Parameters
			----------
			n1 : float
				Refractive index at the wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			R : float
				Atmospheric refraction in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H

		_a = n1 - 1
		_b = H / self.rc
		_R = _a * ( (2 - _a) / (np.sqrt( 2*_b - _a )) ) * np.sin( np.deg2rad(zenith) ) * self.psi( (np.cos(np.deg2rad(zenith))) / np.sqrt(2*_b - _a) )
		return np.degrees(_R)


	def gammaCoefficients(self, n1, reduced_height=None):
		'''
			Submethod for matharExponential(). Calculates the gamma
			coefficients.

			Parameters
			----------
			n1 	: 	float
				Refractive index at the wavelength of interest at the observer
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			g 	: 	array
				List of values, corresponding to the gamma coefficients.
				g[0] corresponds to g_1, g[1] to g_3, g[2] to g_5, and so on.
		'''

		if reduced_height:
			H = reduced_height
		else:
			H = self.H

		_chi0   = n1**2 - 1
		_Khat 	= H / self.rc

		_g1 = (1 - 5/24 * _Khat * _chi0**2 - 3/4 * _chi0 - 3/8 * _Khat**2 * _chi0 \
			- _Khat - 6 * _Khat**3 + 5/8 * _chi0**2 + 63/128 * _chi0**4 + 3/8 \
			* _Khat * _chi0 + 35/256 * _Khat * _chi0**3 + 2 * _Khat**2 + 5/36 \
			* _Khat**2 * _chi0**2 - 35/64 * _chi0**3 + 9/16 * _Khat**3 * _chi0 \
			- 120 * _Khat**5 + 24 * _Khat**4 - 231 / 512 * _chi0**5 - 63 / 640 \
			* _Khat * _chi0**4 - 35 / 512 * _Khat**2 * _chi0**3 - 5 / 36 * _Khat**3\
			* _chi0**2 - 9 / 8 * _Khat**4 * _chi0 )
		
		_g3 = (35/96 * _chi0**3 - 5/12 * _chi0**2 + 333 * _Khat**4 * _chi0 - 2400 \
			* _Khat**5 + 3605/4608 * _Khat**2 * _chi0**3 + 63/256 * _Khat * _chi0**4\
			 - 819/16 * _Khat**3 * _chi0 - 35/18 * _Khat**2 * _chi0**2 - 35/96 \
			 * _Khat * _chi0**3 + 69/8 * _Khat**2 * _chi0 + 5/8 * _Khat * _chi0**2 \
			 - 3/2 * _Khat * _chi0 + 1355/216 * _Khat**3 * _chi0**2 + 10 * _Khat**2 \
			 + 336 * _Khat**4 - 2 * _Khat + _chi0/2 - 21/64 * _chi0**4 - 54 * _Khat**3 \
			 + 77/256 * _chi0**5 ) * 1/2
		
		_g5 = (-96 * _Khat**3 + 8 * _Khat**2 - 7/24 * _chi0**3 + 65/4 * _Khat**2 \
			* _chi0**2 + 35/48 * _Khat * _chi0**3 + 30 * _Khat**2 * _chi0 - 5/3 \
			* _Khat * _chi0**2 - 3 * _Khat * _chi0 + 20091/8 * _Khat**4 * _chi0 \
			- 32605/216 * _Khat**3 * _chi0**2 - 35/9 * _Khat**2 * _chi0**3 + 984 \
			* _Khat**4 + _chi0**2/3 - 10200 * _Khat**5 - 77/320 * _chi0**5 + 21/80 \
			* _chi0**4 - 1089/4 * _Khat**3 * _chi0 - 7/16 * _Khat * _chi0**4 ) * 3/8
		
		_g7 = (-17520 * _Khat**5 - 9/40 * _chi0**4 - 8755/12 * _Khat**3 * _chi0**2 \
			+ 1056 * _Khat**4 + _chi0**3/4 + 23517/4 * _Khat**4 * _chi0 - 11/3 \
			* _Khat * _chi0**2 + 21 * _Khat**2 * _chi0 - 7/4 * _Khat * _chi0**3 \
			+ 505/9 * _Khat**2 * _chi0**2 - 801/2 * _Khat**3 * _chi0 + 33/160 \
			* _chi0**5 - 48 * _Khat**3 + 63/80 * _Khat * _chi0**4 + 2345/96 \
			* _Khat**2 * _chi0**3) * 5/16
		
		_g9 = (_chi0**4 / 5 + 384 * _Khat**4 - 13440 * _Khat**5 + 5562 * _Khat**4 \
			* _chi0 + 3125/36 * _Khat**2 * _chi0**3 - 11/60 * _chi0**5 - 25/6 \
			* _Khat * _chi0**3 - 26950/27 * _Khat**3 * _chi0**2 + 340/9 * _Khat**2 \
			* _chi0**2 - 180 * _Khat**3 * _chi0 - 9/5 * _Khat * _chi0**4) * 35/128

		g = np.array([_g1, _g3, _g5, _g7, _g9])
		g*= n1*_chi0/2

		return g		
	

	def H_isotherm(self, T):
		'''
			Calculates the reduced height of the atmosphere, 
			assuming an isotherm distribution.

			Parameters
			----------
			T 	: float
				Atmospheric temperature in Kelvin

			Returns
			-------
			H 	: float
				Reduced height of the atmosphere in meters
		'''
		_kb = 1.38064852e-23 # kg m2 s-2 K-1
		_m  = 4.809651698e-26 # kg (weight per molecule of air, assuming dry air)

		_phi = np.deg2rad(self.lat)
		_c1 = 5.2790414e-3
		_c2 = 2.32718e-5
		_c3 = 1.262e-7
		_c4 = 7e-10
		_g0 = 9.780327 # m/s^2
		_g0_local = _g0 * (1 + _c1 * np.sin(_phi)**2 + _c2 * np.sin(_phi)**4 + _c3 * np.sin(_phi)**6 + _c4 * np.sin(_phi)**8)

		return (_kb*T) / (_m * _g0_local)


	def matharExponential(self, n1, zenith, reduced_height=None):
		'''
			Model on atmospheric refraction described by Richard Mathar in 
			"A Barometric Exponential Model of the Atmosphere’s Refractive Index:  ZenithAngles and Second Order Aberration in the Entrance Pupil" (2020)
			ArXiv id: 2004.11808 (https://arxiv.org/abs/2004.11808)
			
			The code in this function is adapted from eqs 3, 16, 17 and 18 in this
			work. It is an adaptation of the RefrExpo.cxx C++ file written by Mathar.

			This model assumes the refractive index, n, is related to the susceptibility, 
			chi, following: n(r) = sqrt(1+chi(r)) (eq. 9).

			R. Mathar then assumes that the susceptibility is proportional to the
			air density and reduced height, following eq. 10:
			chi(r) = chi0 * exp(-(r-rho)/K).

			The refraction is described using a Taylor expansion; R in terms of
			the odd powers of the tangent of the observed zenith angle (eq. 3), not
			unlike Oriani's theorem. This requires coefficients g_i.

			Parameters
			----------
			n1 : float
				Refractive index at the wavelength of interest at the observer
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			R : float
				Atmospheric refraction in degrees
		'''

		if reduced_height:
			_H = reduced_height
		else:
			_H = self.H

		_zenith	= np.deg2rad(zenith)
		_rc 	= self.rc # rho

		# Atmospheric refraction in radians
		_gamma  = self.gammaCoefficients(n1, reduced_height=_H)
		_R = sum([_gamma[i]*np.tan(_zenith)**(2*i+1) for i in range(len(_gamma))])
		
		return np.degrees(_R)	


	def oriani(self, n1, zenith, reduced_height=None):
		'''
			Classic refraction formula with two tan(z) terms. This does not 
			assume any information about the structure of the atmosphere. 
			Found with Oriani's theorem.

			Parameters
			----------
			n1 : float
				Refractive index at the wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			R : float
				Atmospheric refraction in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H

		_a = n1 - 1
		_b = H / self.rc
		_R = _a * (1 - _b) * np.tan( np.deg2rad(zenith) ) - _a * (_b - _a/2) * np.tan( np.deg2rad(zenith) )**3
		return np.degrees(_R)


	def planeParallel(self, n1, zenith):
		'''
			Refraction based on a flat atmosphere, easily derived from Snell's law.

			Parameters
			----------
			n1 : float
				Refractive index at the wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.

			Returns
			-------
			R : float
				Atmospheric refraction in degrees
		'''
		_R = np.arcsin(n1 * np.sin( np.deg2rad(zenith) )) - np.deg2rad(zenith)
		return np.degrees(_R)


	def psi(self, x):
		'''
			Submethod for corbard(). Described explicitly in Corbard et al.,2019.

			Parameters
			----------
			x 	: 	float
				Input of the function

			Returns
			-------
			y 	: 	float
				Output of the function
		'''
		_f   = lambda a: np.exp(-a**2)
		_int = integrate.quad(_f, x, np.inf)
		return np.exp(x**2) * _int[0]


	def refractionIntegral(self, n1, zenith, R0=None, useStandardAtmosphere=True,
		heightData=None, rhoData=None):
		'''
			Using data that gives temperature, pressure and density as a 
			function of height, the gladstone-dale relation is invoked to
			calculate the path of a light ray along the atmosphere. 
			This is allows us to also include atmospheric activity in 
			different atmosphere layers.
			Recommended use with the US Standard Atmosphere from 1976.
			USSA1976 data generated from code supplied by http://www.pdas.com/atmos.html
			The method described in eq. 3 in Auer & Standish (2000) is used for the 
			calculation.

			Parameters
			----------
			n1 : float
				Refractive index at the wavelength of interest at the observer
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			R0 : float (optional)
				Local radius of curvature of the Earth in km.
			useStandardAtmosphere : bool
				Flag to use the included USSA1976 data
			heightData : 1D Array
				If "useStandardAtmosphere" is False, the function will ask
				for alternative atmospheric data to integrate over. This
				gives the altitudinal points
			rhoData : 1D Array
				If useStandardAtmosphere is False, the function will ask
				for alternative atmospheric data to integrate over. This
				gives the atmospheric density at the points defined by 
				"heightData".

			Returns
			-------
			R : float
				Atmospheric refraction in degrees

		'''
		if not R0:
			R0 = self.rc
		if useStandardAtmosphere:
			_datafileLocation = os.path.join(os.path.dirname(__file__), 'data/1976USSA.txt')
			_h, _, _, _, _Temp, _Press, _Dens,_,_,_ = np.genfromtxt(_datafileLocation, unpack=True, skip_header=6)
		else:
			_Dens = rhoData
			_h 	  = heightData

		zenith 		= np.deg2rad(zenith)
		_Dens 	= np.asarray(_Dens)
		_h 	 	= np.asarray(_h) 
		_R 		= R0 + _h * 1000
		
		# Gladstone-Dale relation between refractive index and density
		_n_h = 1 + (n1 - 1) / _Dens[0] * _Dens 

		# Calculation of the zenith angle for given height and refractive index, from the refractive invariant
		_z 		= np.arcsin( (n1 * R0) / (_n_h * _R) * np.sin(zenith) )
		
		# Calculation of log space derivative of n and r
		_lnn = np.log(_n_h)
		_lnr = np.log(_R)
		_dlnn = np.asarray([_lnn[i+1] - _lnn[i] for i in range(len(_lnn)-1)])
		_dlnr = np.asarray([_lnr[i+1] - _lnr[i] for i in range(len(_lnr)-1)])
		_dz   = np.asarray([(_z[i] - _z[i+1]) for i in range(len(_z)-1)])

		# Calculation of the integrand in eq. 3 in Auer & Standish
		_Integrand = -1*(_dlnn/_dlnr) / (1 + _dlnn/_dlnr) * _dz


		return np.degrees(np.sum(_Integrand))


	def setReducedHeight(self, p, rho):
		'''
			Sets the reduced atmospheric height.

			Parameters
			----------
			p 	: float
				Atmospheric pressure in Pascal
			rho : float
				Density of the atmospheric air at the observer
		'''
		self.H  	= self.set_H(p, rho) 		# Reduced height of the atmosphere assuming ideal gas law
	

	def set_H(self, p, rho):
		'''
			Calculates the reduced height of the atmosphere, assuming ideal gas law.

			Parameters
			----------
			p 	: float
				Atmospheric pressure in Pascal
			rho : float
				Density of the atmospheric air at the observer

			Returns
			-------
			H 	: float
				Reduced height of the atmosphere in meters
		'''
		_phi = np.deg2rad(self.lat)
		_c1 = 5.2790414e-3
		_c2 = 2.32718e-5
		_c3 = 1.262e-7
		_c4 = 7e-10
		_g0 = 9.780327 # m/s^2
		_g0_local = _g0 * (1 + _c1 * np.sin(_phi)**2 + _c2 * np.sin(_phi)**4 + _c3 * np.sin(_phi)**6 + _c4 * np.sin(_phi)**8)
		_g = _g0_local - (3.0877e-6 - 4.3e-9 * np.sin(_phi)**2) * self.h + 7.2e-13 * self.h**2
		return p / (rho * _g)


	def set_rc(self):
		'''
			Returns the radius of the curvature at the observer in km,
			taking into account the latitude of the observer and the corresponding
			ellipsoidality of the earth.

			Returns
			-------
			rc 	: float
				Local radius of curvature of the Earth in kilometers.
		'''
		_A   = 6378.137 # km
		_B   = 6356.752 # km
		_phi = np.deg2rad(self.lat)
		
		_rc0 = (_A*_B)**2 / (_A**2 * np.cos(_phi)**2 + _B**2 * np.sin(_phi)**2)**(3/2)
		return _rc0 + self.h/1000


	def tan5(self, n1, zenith, reduced_height=None):
		'''
			Same as oriani's theorem, but including a higher order term.
			Found in Corbard et al., 2019.

			Parameters
			----------
			n1 : float
				Refractive index at the wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			R : float
				Atmospheric refraction in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H
		
		_a = n1 - 1
		_b = H / self.rc
		_R = _a * (1 - _b) * np.tan( np.deg2rad(zenith) )  - _a * (_b - _a/2) * np.tan( np.deg2rad(zenith) )**3 + 3 * _a * (_b - _a/2)**2 * np.tan( np.deg2rad(zenith) )**5
		return np.degrees(_R)

