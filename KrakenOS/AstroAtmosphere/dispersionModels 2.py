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
				index of refraction for a wide range of conditions, with corresponding
				error propagation if needed.
				The actual refraction is calculated by any of the following models:
				- Plane parallel atmosphere
				- Cassini spherical homogeneous atmosphere
				- Mathar's Barometric exponential model
				- Other methods described in Corbard et al., 2019.
'''

class dispersion(Observatory):
	'''
		Class that sets the observer location.
		These location attributes are then used to calculate the 
		appropriate atmospheric dispersion.

		See also:
		Observatory()
		refraction()
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
			Returns the atmospheric dispersion assuming Cassini's homegeneous
			atmosphere model.
		cassiniError(z, l1, l2, T, p, RH, xc, dl1=0, dl2=0, dT=0.2, dP=20, 
					dRH=0.02, dCO2=20, dz=0, lat=None, h=None)
			Returns the uncertainty in the atmospheric dispersion, assuming 
			Cassini's homegeneous atmosphere model and the given parameters.
		corbard(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric dispersion following the method described 
			in Corbard et al., 2019.
		gammaCoefficients(n1, reduced_height=None)
			Helper function for matharExponential(), that calculates the gamma
			coefficients.
		H_isotherm(T)
			Sets the reduced height for the conditions at the observer, assuming
			an isothermal atmosphere.
		matharExponential(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric dispersion assuming	Mathar's barometric 
			exponential model.	
		oriani(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric dispersion following the classic version
			of Oriani's theorem.
		planeParallel(n1, n2, zenith)
			Returns the atmospheric dispersion assuming a plane parallel
			atmosphere.	
		psi(x)
			Submethod for corbard().
		refractionIntegral(n1, n2, zenith, R0=None, useStandardAtmosphere=True, heightData=None, rhoData=None)
			Returns the atmospheric dispersion following the refraction integral.	
		setReducedHeight(p, rho)
			Sets the reduced atmospheric height.
		set_H(p, rho)
			Sets the reduced height for the conditions at the observer, assuming
			an adiabatic atmosphere.
		set_rc()
			Sets the local radius of curvature of the Earth at the observer.
		tan5(n1, n2, zenith, reduced_height=None)
			Returns the atmospheric dispersion following Oriani's theorem,
			including a fifth order term in the expansion.		
		Other Observatory() methods
	'''
	def __init__(self, lat, h):
		Observatory.__init__(self)
		self.h  	= h 						# height above sea level in meters
		self.lat 	= lat						# Latitude of the observer in degrees
		self.rc 	= self.set_rc() * 1000		# Radius of curvature of the earth at the observer


	def cassini(self, n1, n2, zenith, reduced_height=None):
		'''
			Refraction of spherical atmosphere, derived from geometry using Sine law

			Parameters
			----------
			n1 : float
				Refractive index at the shortest wavelength of interest
			n2 : float
				Refractive index at the longest wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			dR : float
				Atmospheric dispersion in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H
		
		_R1 = np.arcsin(n1 * self.rc * np.sin(np.deg2rad(zenith)) / (self.rc + H)) - np.arcsin(self.rc * np.sin(np.deg2rad(zenith)) / (self.rc + H) )
		_R2 = np.arcsin(n2 * self.rc * np.sin(np.deg2rad(zenith)) / (self.rc + H)) - np.arcsin(self.rc * np.sin(np.deg2rad(zenith)) / (self.rc + H) )	
		return np.degrees(_R1) - np.degrees(_R2)


	def cassiniError(self, zenith, l1, l2, T, p, RH, xc, dl1=0, dl2=0, dT=0.2, dP=20,
						dRH=0.02, dCO2=20, dz=0, lat=None, h=None):
		'''
			Calculates the uncertainty in the atmospheric dispersion, 
			following the Cassini model, by using fully propagated error analysis.
			This involves almost endless use of the chain rule and evaluation 
			of both the moist air and dry air components. 
			Note: dry air does not contain any water, i.e. H=0, x_w=0, etc.

			Not taken into account is the errror propagation in the scale height
			of the atmosphere, from my results this should give a small or even 
			negligible difference.

			Parameters
			----------
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			l1 	: float
				Wavelength of light in microns, at the shortest wavelength 
				of interest
			l2 	: float
				Wavelength of light in microns, at the longest wavelength 
				of interest
			T 	: float
				Atmospheric temperature in Kelvin
			p 	: float
				Atmospheric pressure in Pascal
			RH 	: float
				Atmospheric relative humidity (0<=RH<=1)
			xc 	: float
				CO2 density in parts per million
			dl1	: float 
				Uncertainty in the wavelength in microns, at the shortest
				wavelength of interest
			dl2	: float 
				Uncertainty in the wavelength in microns, at the longest
				wavelength of interest
			dT 	: float (default: 0.2)
				Uncertainty in the atmospheric temperature in Kelvin
			dP 	: float (default: 20)
				Uncertainty in the atmospheric pressure in Pascal
			dRH : float (default: 0.02)
				Uncertainty in the atmospheric relative humidity
			dCO2: float (default: 20)
				Uncertainty in the atmospheric CO2 density in parts per million
			dz 	: float (default: 0)
				Uncertainty in the zenith angle, in degrees.
			lat : float (default: None)
				If entered, overrides the latitude attribute of the class.
			h 	: float (default=None)
				If entered, overrides the altitude attribute of the class.

			Returns
			-------
			dRR	: float
				Uncertainty in the atmospheric dispersion in degrees.
		'''
		if lat == None:
			lat = self.lat
		else:
			self.lat = lat
		if h == None:
			h 	= self.h
		else:
			self.h = h

		# constants
		_cf = 1.022
		_w0 = 295.235 	# um^-2
		_w1 = 2.6422	# um^-2
		_w2 = -0.032380	# um^-4
		_w3 = 0.004028	# um^-6

		_k0 = 238.0185 		# um^-2
		_k1 = 5792105  		# um^-2
		_k2 = 57.362		# um^-2
		_k3 = 167917		# um^-2	

		_A =  1.2378847e-5 # K^-2
		_B = -1.9121316e-2 # K^-2
		_C = 33.93711047
		_D = -6.3431645e3  # K^-2

		_Mw  = 0.018015   # kg/mol
		_R   = 8.314510	  # J mol^-1 K^-1, gas constant

		_a0 =  1.58123e-6 # K Pa^-1
		_a1 = -2.9331e-8  # Pa^-1
		_a2 =  1.1043e-10 # K^-1 Pa^-1
		_b0 =  5.707e-6   # K Pa^-1
		_b1 = -2.051e-8   # Pa^-1
		_c0 =  1.9898e-4  # K Pa^-1
		_c1 = -2.376e-6   # Pa^-1
		_d  =  1.83e-11   # K^2 Pa^-2
		_e  = -0.765e-8   # K^2 Pa^-2
		_t  = T - 273.15

		_aa = 1.00062
		_bb = 3.14e-8 	# Pa^-1
		_cc = 5.6e-7	# deg C^-2

		_RH_a = 0
		_RH_w = RH

		# Calculated values
		_rho_obs = self.rho(p, T, RH, xc)
		self.rc  = self.set_rc() * 1000
		self.H   = self.set_H(p, _rho_obs) 
		_alpha   = self.rc * np.sin(np.deg2rad(zenith)) / (self.rc + self.H)
		_p_a 	 = p - RH*self.svp(T) 
		_p_w     = RH*self.svp(T)
		_x_w_a	 = self.xw(_p_a, T, _RH_a)
		_x_w_w 	 = self.f(_p_w, T)
		_Z_a 	 = self.Z(p=_p_a, T=T, RH=_RH_a)
		_Z_w 	 = self.Z(p=_p_w, T=T, RH=_RH_w, xw=_x_w_w)
		_Ma  	 = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400) )	# Molar mass in kg/mol of dry air containing xc ppm of CO2
		_p_sv 	 = self.svp(T)
		_rho_a 	 = self.rho(p=p-RH*self.svp(T), T=T, RH=0, xc=xc)
		_rho_axs = self.rho(p=101325, T=288.15, RH=0, xc=xc)
		if RH == 0.0:
			# nan is returned for _rho_w if we put in a (partial) pressure of 0
			_rho_w = 0
		else:
			_rho_w = self.rho(p=RH*self.svp(T), T=T, RH=RH, xc=xc) 
		_rho_ws  = self.rho(p=1333, T=293.15, RH=1333/self.svp(293.15), xc=xc)
		_n_a1 	 = self.n_axs(l1, xc)
		_n_w1 	 = self.n_ws(l1)
		_n_a2 	 = self.n_axs(l2, xc)
		_n_w2 	 = self.n_ws(l2)

		_L1_a 	 = (_n_a1**2 - 1) / (_n_a1**2 + 2)  # Dry air component of L
		_L1_w 	 = (_n_w1**2 - 1) / (_n_w1**2 + 2)
		_L1   	 = ( _rho_a / _rho_axs ) * _L1_a + ( _rho_w / _rho_ws ) * _L1_w # Overall L
		_L2_a 	 = (_n_a2**2 - 1) / (_n_a2**2 + 2) # Dry air component of L
		_L2_w 	 = (_n_w2**2 - 1) / (_n_w2**2 + 2)
		_L2    	 = ( _rho_a / _rho_axs ) * _L2_a + ( _rho_w / _rho_ws ) * _L2_w # Overall L

		_n1   	 = np.sqrt( (1 + 2 * _L1) / (1 - _L1) )
		_n2   	 = np.sqrt( (1 + 2 * _L2) / (1 - _L2) )

		# Partial derivatives for derivation to l1
		_dn_wdl1 = _cf * 1e-8 * (-1 * _w1 * 2/l1**3 - _w2 * 4/l1**5 - _w3 * 6/l1**7)
		_dLdn_w1 = _rho_w / _rho_ws * 6 * _n_w1 / (_n_w1**2 + 2)**2
		_dn_adl1 = 1e-8 * (1 + 0.534e-6 * (xc - 450)) * ( - (2 * _k1 * l1 ) / ((_k0 * l1**2 - 1)**2) - (2 * _k3 * l1) / ((_k2 * l1**2 - 1)**2) )
		_dLdn_a1 = _rho_a / _rho_axs * 6 * _n_a1 / (_n_a1**2 + 2)**2
		_dndL1   = ( ( (1 + 2*_L1) / (1 - _L1)**2 ) + (2 / (1-_L1)) ) / (2 * np.sqrt( (1 + 2*_L1) / (1 - _L1) ) ) 
		_dddn1   = _alpha / np.sqrt( 1 - _alpha**2 * _n1**2 )
		_dddl1   = _dddn1 * _dndL1 * ( _dLdn_a1 * _dn_adl1 + _dLdn_w1 * _dn_wdl1)

		# Partial derivatives for derivation to l2
		_dn_wdl2 = _cf * 1e-8 * (-1 * _w1 * 2/l2**3 - _w2 * 4/l2**5 - _w3 * 6/l2**7)
		_dLdn_w2 = _rho_w / _rho_ws * 6 * _n_w2 / (_n_w2**2 + 2)**2
		_dn_adl2 = 1e-8 * (1 + 0.534e-6 * (xc - 450)) * ( - (2 * _k1 * l2 ) / ((_k0 * l2**2 - 1)**2) - (2 * _k3 * l2) / ((_k2 * l2**2 - 1)**2) )
		_dLdn_a2 = _rho_a / _rho_axs * 6 * _n_a2 / (_n_a2**2 + 2)**2
		_dndL2   = ( ( (1 + 2*_L2) / (1 - _L2)**2 ) + (2 / (1-_L2)) ) / (2 * np.sqrt( (1 + 2*_L2) / (1 - _L2) ) ) 
		_dddn2   = _alpha / np.sqrt( 1 - _alpha**2 * _n2**2 )
		_dddl2   = _dddn2 * _dndL2 * ( _dLdn_a2 * _dn_adl2 + _dLdn_w2 * _dn_wdl2)


		# Partial derivatives for derivation to T
		# Dry air parts:
		_drho_adp_a  = _Ma / (_Z_a * _R * T) * (1 - _x_w_a * (1 - _Mw / _Ma))
		_dp_adp_sv 	 = - RH
		_dp_svdT 	 = (2*_A*T + _B - _D/(T**2)) * np.exp(_A * T**2 + _B * T + _C + _D/T )
		_drho_adTZ	 = - _p_a * _Ma / (_R * ( _Z_a * T)**2) * (1 - _x_w_a * (1 - _Mw / _Ma))
		_dTZdT_a	 = _p_a / T * (_a2 * _t**2 + _a1 * _t + _a0 + (_b1*_t + _b0) * _x_w_a + (_c0 + _c1*_t) * _x_w_a**2 ) - _p_a * (_a1 + 2*_a2*_t + _b1*_x_w_a + _c1*_x_w_a**2) - 2*(_p_a/T)**2 * (_d + _e * _x_w_a**2) + _Z_a
		_dTZdxw_a	 = 2*_e*_p_a**2 / T * _x_w_a - _p_a * (_b1*_t + _b0 + 2 * (_c0 + _c1*_t) * _x_w_a)
		_dxw_adf	 = _RH_a * _p_sv / _p_a
		_dfdT 		 = 2*_cc*_t
		_dxw_adp_sv	 = self.f(_p_a, T) * _RH_a / _p_a
		_drho_adxw_a = 0 # as rho_a has not water content, i.e. rho = p*M/ZRT


		# Moist air parts:
		_drho_wdp_w  = _Ma / (_Z_w * _R * T) * (1 - _x_w_w * (1 - _Mw / _Ma))
		_dp_wdp_sv 	 = RH
		_dp_svdT 	 = (2*_A*T + _B - _D/(T**2)) * np.exp(_A * T**2 + _B * T + _C + _D/T )
		_drho_wdTZ	 = - _p_w * _Ma / (_R * (_Z_w * T)**2) * (1 - _x_w_w * (1 - _Mw / _Ma))
		_dTZdT_w	 = _p_w / T * (_a2 * _t**2 + _a1 * _t + _a0 + (_b1*_t + _b0) * _x_w_w + (_c0 + _c1*_t) * _x_w_w**2 ) - _p_w * (_a1 + 2*_a2*_t + _b1*_x_w_w + _c1*_x_w_w**2) - 2*(_p_w/T)**2 * (_d + _e * _x_w_w**2) + _Z_w
		_dTZdxw_w	 = 2*_e*_p_w**2 / T * _x_w_w - _p_w * (_b1*_t + _b0 + 2 * (_c0 + _c1*_t) * _x_w_w)
		_dxw_wdf	 = 1 # _RH_w * _p_sv / _p_w
		_dfdT 		 = 2*_cc*_t
		_dxw_wdp_sv	 = self.f(_p_w, T) / _p_sv # * _RH_w / _p_w
		_drho_wdxw_w = _p_w * _Ma / (_Z_w * _R * T) * (_Mw / _Ma - 1)

		_tA =  _drho_adp_a * _dp_adp_sv * _dp_svdT +  _drho_adTZ * (_dTZdT_a + _dTZdxw_a * (_dxw_adf * _dfdT + _dxw_adp_sv * _dp_svdT) ) + _drho_adxw_a * ( _dxw_adf * _dfdT + _dxw_adp_sv * _dp_svdT )
		_tB =  _drho_wdp_w * _dp_wdp_sv * _dp_svdT +  _drho_wdTZ * (_dTZdT_w + _dTZdxw_w * (_dxw_wdf * _dfdT + _dxw_wdp_sv * _dp_svdT) ) + _drho_wdxw_w * ( _dxw_wdf * _dfdT + _dxw_wdp_sv * _dp_svdT )
		_dL1drho_a 	 = _L1_a / _rho_axs
		_dL2drho_a 	 = _L2_a / _rho_axs
		_dL1drho_w 	 = _L1_w / _rho_ws
		_dL2drho_w 	 = _L2_w / _rho_ws
		
		_dddT = _tA * (_dddn1 * _dndL1 * _dL1drho_a - _dddn2 * _dndL2 * _dL2drho_a) + _tB * (_dddn1 * _dndL1 * _dL1drho_w - _dddn2 * _dndL2 * _dL2drho_w)

		# Partial derivatives for derivation to p
		# Dry air parts:
		_drho_adZ 	= - _p_a * _Ma / ( _Z_a**2 * _R * T) * (1 - _x_w_a * (1 - _Mw / _Ma))
		_dZdp_a 	= - 1 / T * (_a0 + _a1 * _t + _a2 * _t**2 + (_b0 + _b1*_t)*_x_w_a + (_c0 + _c1*_t)*_x_w_a**2 ) + 2*_p_a/(T**2) * (_d + _e * _x_w_a**2)
		_dfdp 		= _bb
		_dxw_adp_a 	= - self.f(_p_a, T) * _RH_a * _p_sv / (_p_a)**2 

		# Moist air parts
		_dZdp_w 	= - 1 / T * (_a0 + _a1 * _t + _a2 * _t**2 + (_b0 + _b1*_t)*_x_w_w + (_c0 + _c1*_t)*_x_w_w**2 ) + 2*_p_w/(T**2) * (_d + _e * _x_w_w**2)
		# Not necessary, since dpwdp=0

		_pA = _drho_adp_a + _drho_adZ * _dZdp_a + _drho_adxw_a * (_dxw_adf * _dfdp + _dxw_adp_a)
		_dddP = _pA * (_dddn1 * _dndL1 * _dL1drho_a - _dddn2 * _dndL2 * _dL2drho_a)

		
		# Partial derivatives for derivation to H
		# Dry air parts:
		_dp_adH 	= - _p_sv
		_dZdxw_a 	= (_p_a / T)**2 * 2 * _e * _x_w_a - (_p_a / T) * ( (_b0 + _b1 * _t) + 2 * (_c0 + _c1 * _t) * _x_w_a )
		_dxw_adH 	= self.f(_p_a, T) * _p_sv / _p_a

		# Moist air parts:
		_dp_wdH 	= _p_sv
		_dZdxw_w 	= (_p_w / T)**2 * 2 * _e * _x_w_w - (_p_w / T) * ( (_b0 + _b1 * _t) + 2 * (_c0 + _c1 * _t) * _x_w_w )
		_dxw_wdH 	= 0
		_dxw_wdp_w 	= 0		
		_drho_wdZ 	= - _p_w * _Ma / ( _Z_w**2 * _R * T) * (1 - _x_w_w * (1 - _Mw / _Ma))

		_hA = _drho_adp_a * _dp_adH + _drho_adZ * ( _dZdp_a * _dp_adH + _dZdxw_a * (_dxw_adH + _dxw_adp_a * _dp_adH + _dxw_adf * _dfdp * _dp_adH) ) + _drho_adxw_a * (_dxw_adH + _dxw_adp_a * _dp_adH + _dxw_adf * _dfdp * _dp_adH)
		_hB = _drho_wdp_w * _dp_wdH + _drho_wdZ * ( _dZdp_w * _dp_wdH + _dZdxw_w * ( _bb * _p_sv ) ) + _drho_wdxw_w * ( _bb * _p_sv )
		_dddRH = _hA * (_dddn1 * _dndL1 * _dL1drho_a - _dddn2 * _dndL2 * _dL2drho_a) + _hB * (_dddn1 * _dndL1 * _dL1drho_w - _dddn2 * _dndL2 * _dL2drho_w)


		# Partial derivatives for derivation to xc
		# Dry air parts:
		_drho_adMa   = _p_a * (1 - _x_w_a) / (_Z_a * _R * T)
		_dL1drho_axs = - (_rho_a * _L1_a)/(_rho_axs**2)
		_dL2drho_axs = - (_rho_a * _L2_a)/(_rho_axs**2)
		_drho_axsdMa = _drho_adMa
		_dL1dL_a	 = _rho_a / _rho_axs
		_dL_adn_a1	 = 6 * _n_a1 / (_n_a1**2 + 2)**2
		_dL2dL_a	 = _dL1dL_a
		_dL_adn_a2	 = 6 * _n_a2 / (_n_a2**2 + 2)**2

		# Moist air parts:
		_drho_wdMa   = _p_w * (1 - _x_w_w) / (_Z_w * _R * T)
		_dL1drho_ws = - (_rho_w * _L1_w)/(_rho_ws**2)
		_dL2drho_ws = - (_rho_w * _L2_w)/(_rho_ws**2)
		_drho_wsdMa = _drho_wdMa

		# Other parts:
		_dMadx_c   = 12.011e-9
		_dn_adx_c1 = (self.n_as(l1) - 1) * 0.534e-6
		_dn_adx_c2 = (self.n_as(l2) - 1) * 0.534e-6

		_dL1dx_c = _dL1drho_a * _drho_adMa * _dMadx_c + _dL1drho_w * _drho_wdMa * _dMadx_c + _dL1drho_axs * _drho_axsdMa * _dMadx_c + _dL1drho_ws * _drho_wsdMa * _dMadx_c + _dL1dL_a * _dL_adn_a1 * _dn_adx_c1
		_dL2dx_c = _dL2drho_a * _drho_adMa * _dMadx_c + _dL2drho_w * _drho_wdMa * _dMadx_c + _dL2drho_axs * _drho_axsdMa * _dMadx_c + _dL2drho_ws * _drho_wsdMa * _dMadx_c + _dL2dL_a * _dL_adn_a2 * _dn_adx_c2
		_dddCO2 = _dddn1 * _dndL1 * _dL1dx_c - _dddn2 * _dndL2 * _dL2dx_c


		# Partial derivatives for derivation to z
		_zenith = np.deg2rad(zenith)
		_dz   	= np.deg2rad(dz)
		_R 		= self.rc
		_H 		= self.H
		_dddz 	= (_R * np.cos(_zenith)) / (_R + _H) * ( _n1 / np.sqrt(1-(_n1**2 * _R**2 * np.sin(_zenith)**2)/(_R + _H)**2) - _n2 / np.sqrt(1-(_n2**2 * _R**2 * np.sin(_zenith)**2)/(_R + _H)**2) )		

		# Final error:
		_dispersionError = np.sqrt( (_dddl1 * dl1)**2 + (_dddl2 * dl2)**2 + (_dddT * dT)**2  + (_dddP * dP)**2 + (_dddRH * dRH)**2 + (_dddCO2 * dCO2)**2  + (_dddz * _dz)**2)

		return np.degrees(_dispersionError)


	def corbard(self, n1, n2, zenith, reduced_height=None):
		'''
			Corbard et al., 2019, mention an additional formula based on
			the error function. Oriani's theorem can be derived from this
			equation by 'keeping only the three first terms of its 
			asymptotic expansion'.

			Parameters
			----------
			n1 : float
				Refractive index at the shortest wavelength of interest
			n2 : float
				Refractive index at the longest wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			dR : float
				Atmospheric dispersion in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H

		_rList = []
		for _n in [n1, n2]:
			_a = _n - 1
			_b = H / self.rc
			_R = _a * ( (2 - _a) / (np.sqrt( 2*_b - _a )) ) * np.sin( np.deg2rad(zenith) ) * self.psi( (np.cos(np.deg2rad(zenith))) / np.sqrt(2*_b - _a) )
			_rList.append(np.degrees(_R))
		return _rList[0] - _rList[1]


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


	def matharExponential(self, n1, n2, zenith, reduced_height=None):
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
			n2 : float
				Refractive index at the longest wavelength of interest at the
				observer
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			dR : float
				Atmospheric refraction in degrees
		'''
		if reduced_height:
			_H = reduced_height
		else:
			_H = self.H

		_zenith	= np.deg2rad(zenith)
		_rc 	= self.rc # rho

		# Atmospheric refraction for the first wavelength
		_gamma1 = self.gammaCoefficients(n1, reduced_height=_H)
		_R1 = sum([_gamma1[i]*np.tan(_zenith)**(2*i+1) for i in range(len(_gamma1))])

		# Atmospheric refraction for the second wavelength
		_gamma2 = self.gammaCoefficients(n2, reduced_height=_H)
		_R2 = sum([_gamma2[i]*np.tan(_zenith)**(2*i+1) for i in range(len(_gamma2))])

		# The atmospheric dispersion in radians
		_dR = _R1 - _R2

		return np.degrees(_dR)


	def oriani(self, n1, n2, zenith, reduced_height=None):
		'''
			Classic refraction formula with two tan(z) terms. 
			This does not assume any information about the structure
			of the atmosphere. Found with Oriani's theorem.

			Parameters
			----------
			n1 : float
				Refractive index at the shortest wavelength of interest
			n2 : float
				Refractive index at the longest wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			dR : float
				Atmospheric dispersion in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H

		_rList = []
		for _n in [n1, n2]:
			_a = _n - 1
			_b = H / self.rc
			_R = _a * (1 - _b) * np.tan( np.deg2rad(zenith) ) - _a * (_b - _a/2) * np.tan( np.deg2rad(zenith) )**3
			_rList.append(np.degrees(_R))
		return _rList[0] - _rList[1]


	def planeParallel(self, n1, n2, zenith):
		'''
			Refraction based on a flat atmosphere, easily derived from Snell's law.

			Parameters
			----------
			n1 : float
				Refractive index at the shortest wavelength of interest
			n2 : float
				Refractive index at the longest wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.

			Returns
			-------
			dR : float
				Atmospheric dispersion in degrees
		'''
		_R1 = np.arcsin(n1 * np.sin( np.deg2rad(zenith) )) - np.deg2rad(zenith)
		_R2 = np.arcsin(n2 * np.sin( np.deg2rad(zenith) )) - np.deg2rad(zenith)
		return np.degrees(_R1) - np.degrees(_R2)


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


	def refractionIntegral(self, n1, n2, zenith, R0=None, useStandardAtmosphere=True, 
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
				Refractive index at the shortest wavelength of interest
			n2 : float
				Refractive index at the longest wavelength of interest
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
			dR : float
				Atmospheric dispersion in degrees
		'''
		if not R0:
			R0 = self.rc
		if useStandardAtmosphere:
			_datafileLocation = os.path.join(os.path.dirname(__file__), 'data/1976USSA.txt')
			_h, _, _, _, _Temp, _Press, _Dens,_,_,_ = np.genfromtxt(_datafileLocation, unpack=True, skip_header=6)
		else:
			_Dens = rhoData
			_h 	  = heightData 

		
		_rList 	= []
		zenith 	= np.deg2rad(zenith)
		
		# Integration of the refraction integral over the atmosphere layers.
		for _n in [n1, n2]:
			_Dens 	= np.asarray(_Dens)
			_h 	 	= np.asarray(_h) 
			_R 		= R0 + _h * 1000
			
			# Gladstone-Dale relation between refractive index and density
			_n_h = 1 + (_n - 1) / _Dens[0] * _Dens 

			# Calculation of the zenith angle for given height and refractive index, from the refractive invariant
			_z 		= np.arcsin( (_n * R0) / (_n_h * _R) * np.sin(zenith) )
			
			# Calculation of log space derivative of n and r
			_lnn = np.log(_n_h)
			_lnr = np.log(_R)
			_dlnn = np.asarray([_lnn[i+1] - _lnn[i] for i in range(len(_lnn)-1)])
			_dlnr = np.asarray([_lnr[i+1] - _lnr[i] for i in range(len(_lnr)-1)])
			_dz   = np.asarray([(_z[i] - _z[i+1]) for i in range(len(_z)-1)])

			# Calculation of the integrand in eq. 3 in Auer & Standish
			_Integrand = -1*(_dlnn/_dlnr) / (1 + _dlnn/_dlnr) * _dz
			_rList.append(np.degrees(np.sum(_Integrand)))
		return _rList[0] - _rList[1]


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


	def tan5(self, n1, n2, zenith, reduced_height=None):
		'''
			Same as oriani's theorem, but including a higher order term.
			Found in Corbard et al., 2019.

			Parameters
			----------
			n1 : float
				Refractive index at the shortest wavelength of interest
			n2 : float
				Refractive index at the longest wavelength of interest
			zenith : float
				Angle of observation in degrees.
				zenith = 0 means the observation is directly overhead.
			reduced_height : float (optional)
				Reduced height of the atmosphere. Default value is calculated
				from the object attributes.

			Returns
			-------
			dR : float
				Atmospheric dispersion in degrees
		'''
		if reduced_height:
			H = reduced_height
		else:
			H = self.H
		
		_rList = []
		for _n in [n1, n2]:
			_a = _n - 1
			_b = H / self.rc
			_R = _a * (1 - _b) * np.tan( np.deg2rad(zenith) )  - _a * (_b - _a/2) * np.tan( np.deg2rad(zenith) )**3 + 3 * _a * (_b - _a/2)**2 * np.tan( np.deg2rad(zenith) )**5
			_rList.append(np.degrees(_R))
		return _rList[0] - _rList[1]
