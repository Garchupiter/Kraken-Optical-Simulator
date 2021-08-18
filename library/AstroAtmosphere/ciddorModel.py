#!/usr/bin/env Python3
import numpy as np
from scipy import integrate

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
				- Other methods described in Corbard et al., 2019.
'''

class Observatory:
	'''
		Class used for all function to calculate the refractive index
		and the corresponding uncertainties.

		...

		Attributes
		----------
		dT 	: 	float (optional)
			Uncertainty in the atmospheric temperature in Kelvin
		dP 	: 	float (optional)
			Uncertainty in the atmospheric pressure in Pascal
		dH 	: 	float (optional)
			Uncertainty in the atmospheric relative humidity
		dCO2: 	float (optional)
			Uncertianty in the atmospheric CO2 density in parts
			per million.
		dl  : 	float (optional)
			Uncertainty in the wavelengths in micron

		Methods
		-------
		n_as(l)
			Returns the refractive index for standard atmospheric conditions.
		n_axs(l, xc)
			Same as n_as(l), but with a different CO2 density given by xc
		n_ws(l)
			Returns the refractive index of air for standard moist contions.
		f(p,T)
			Returns the enchancement factor of water vapour in air.
		svp(T)
			Returns the saturation vapour pressure in Pa.
		svp2(T)
			Returns the saturation vapour pressure in Pa. 
			From IAPWS formula, based on documentation by NIST toolbox.
		xw(p, T, RH)
			Returns the molar fraction of water vapour in moist air.
		Z(p, T, RH, xw='')
			Returns the compressibility of moist air
		rho(p, T, RH, xc, xw='')
			Returns the density of atmospheric air at the given parameters.
		n_tph(l, T, p, RH, xc)
			Returns the refractive index for any atmospheric conditions.
		dn_tph(l, T, p, RH, xc, dl=None, dT=None, dP=None, dRH=None, dCO2=None)
			Returns the uncertainty in the refractive index for any atmospheric
			conditions, given the parameters.
	'''
	def __init__(self, dT=0.2, dP=20, dH=0.02, dCO2=20, dl=0.001):
		self.dT 	= dT 	# in K
		self.dP 	= dP 	# in Pa
		self.dH 	= dH	# relative humidity between 0 and 1
		self.dCO2 	= dCO2	# in ppm
		self.dl 	= dl 	# in micron

	def n_as(self, l):
		'''
			Formula that calculate the refractive index for standard atmospheric
			conditions, i.e. 15 C, 101,325 Pa, 0% RH, 450 ppm CO2.
			It returns n_as, not (n_as - 1) * 1e8!

			Parameters
			----------
			l 	: float
				Wavelength of light in microns

			Returns
			-------
			n 	: float
				Refractive index of air at standard atmospheric conditions.

		'''
		_s = 1 / l
		_k0 = 238.0185 	# um^-2
		_k1 = 5792105  	# um^-2
		_k2 = 57.362		# um^-2
		_k3 = 167917		# um^-2
		return (_k1 / (_k0 - _s**2) + _k3 / (_k2 - _s**2)) / 1e8 + 1


	def n_axs(self, l, xc):
		'''
			Formula that calculate the refractive index for standard atmospheric
			conditions, i.e. 15 C, 101,325 Pa, 0% RH, but with xc ppm of CO2.
			It returns n_axs, not (n_axs - 1)

			Parameters
			----------
			l 	: float
				Wavelength of light in microns
			xc 	: float
				CO2 density in parts per million

			Returns
			-------
			n 	: float
				Refractive index of air
		'''
		_n_as = self.n_as(l)
		return ( (_n_as - 1) * (1 + 0.534e-6 * (xc - 450) ) ) + 1


	def n_ws(self, l):
		'''
			Formula that calculates the refractivity of air for standard moist conditions
			as defined by B&S (1939), i.e. 20 C, 1333 Pa
			Returns n_ws, not (n_ws - 1) * 1e8

			Parameters
			----------
			l 	: float
				Wavelength of light in microns

			Returns
			-------
			n 	: float
				Refractive index of standard moist air.
		'''
		_s  = 1 / l
		_cf = 1.022
		_w0 = 295.235 	# um^-2
		_w1 = 2.6422	# um^-2
		_w2 = -0.032380	# um^-4
		_w3 = 0.004028	# um^-6
		return (_cf * (_w0 + _w1 * _s**2 + _w2 * _s**4 + _w3 * _s**6) ) / 1e8 + 1


	def f(self, p, T):
		'''
			Enhancement factor of water vapour in air

			Parameters
			----------
			p 	: float
				Atmospheric pressure in Pascal
			T 	: float
				Atmospheric temperature in Kelvin

			Returns
			-------
			f 	: float
				Enhancement factor of water vapour in air
		'''
		_t = T - 273.15
		_a = 1.00062
		_b = 3.14e-8 	# Pa^-1
		_c = 5.6e-7		# deg C^-2
		return _a + _b * p + _c * _t**2


	def svp(self, T):
		'''
			Saturation vapour pressure in Pa.

			Parameters
			----------
			T 	: float
				Atmospheric temperature in Kelvin

			Returns
			-------
			svp	: float
				Saturation vapour pressure
		'''
		_A =  1.2378847e-5 # K^-2
		_B = -1.9121316e-2 # K^-2
		_C = 33.93711047
		_D = -6.3431645e3  # K^-2
		return np.exp(_A*T**2 + _B*T + _C + _D/T)


	def svp2(self, T):
		'''
			Calculation of svp from IAPWS formula, based on documentation 
			by NIST toolbox.

			Parameters
			----------
			T 	: float
				Atmospheric temperature in Kelvin

			Returns
			-------
			svp	: float
				Saturation vapour pressure
		'''
		K1 =  1.16705214528e+03
		K2 = -7.24213167032e+05
		K3 = -1.70738469401e+01
		K4 =  1.20208247025e+04
		K5 = -3.23255503223e+06
		K6 =  1.49151086135e+01
		K7 = -4.82326573616e+03
		K8 =  4.05113405421e+05
		K9 = -2.38555575678e-01
		K10 = 6.50175348448e+02
		_W = T + K9 / (T - K10)
		_A = _W**2 + K1*_W + K2
		_B = K3*_W**2 + K4*_W + K5
		_C = K6*_W**2 + K7*_W + K8
		_X = -_B + np.sqrt(_B**2 - 4*_A*_C)
		return (2*_C / _X)**4 * 1e6


	def xw(self, p, T, RH):
		'''
			Molar fraction of water vapour in moist air.

			Parameters
			----------
			p 	: float
				Atmospheric pressure in Pascal
			T 	: float
				Atmospheric temperature in Kelvin
			RH 	: float
				Atmospheric relative humidity (0<=RH<=1)

			Returns
			-------
			xw	: float
				Molar fraction of water vapour in moist air
		'''
		_f 		= self.f(p, T)
		_svp 	= self.svp(T)
		return _f * RH * _svp / p


	def Z(self, p, T, RH, xw=''):
		'''
			Compressibility of moist air

			Parameters
			----------
			p 	: float
				Atmospheric pressure in Pascal
			T 	: float
				Atmospheric temperature in Kelvin
			RH 	: float
				Atmospheric relative humidity (0<=RH<=1)

			Returns
			-------
			Z 	: float
				Compressibility of moist air
		'''
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
		if xw == '':
			_xw = self.xw(p, T, RH)
		else:
			_xw = xw
		Z = 1 - (p/T) * (_a0 + _a1*_t + _a2*_t**2 + (_b0 + _b1*_t)*_xw + (_c0 + _c1*_t)*_xw**2) + (p/T)**2 * (_d + _e*_xw**2)
		return Z


	def rho(self, p, T, RH, xc, xw=''):
		'''
			Calculates the density of air for the given parameters.
			
			Parameters
			----------
			p 	: float
				Atmospheric pressure in Pascal
			T 	: float
				Atmospheric temperature in Kelvin
			RH 	: float
				Atmospheric relative humidity (0<=RH<=1)
			xc 	: float
				CO2 density in parts per million
			xw 	: str or float
				Placeholder for molar fraction of water vapour in moist air.
				Used to pass the function along where it should go.

			Returns
			-------
			rho : float
				Density of air
		'''
		_Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400) )	# Molar mass in kg/mol of dry air containing xc ppm of CO2
		_Mw = 0.018015 										# kg/mol
		_R  = 8.314510										# J mol^-1 K^-1, gas constant
		if xw == '':
			_xw = self.xw(p, T, RH)
		else: 
			_xw = xw
		_Z  = self.Z(p, T, RH, xw=_xw)
		return ( ( p * _Ma ) / (_Z * _R * T) ) * ( 1 - _xw * (1 - _Mw/_Ma)) 


	def n_tph(self, l, T, p, RH, xc):
		'''
			Calculates the refraction of the atmosphere for any conditions.
			Follows the equations by Ciddor (1996).

			Parameters
			----------
			l 	: float
				Wavelength of light in microns
			p 	: float
				Atmospheric pressure in Pascal
			T 	: float
				Atmospheric temperature in Kelvin
			RH 	: float
				Atmospheric relative humidity (0<=RH<=1)
			xc 	: float
				CO2 density in parts per million

			Returns
			-------
			n 	: float
				Refractive index of air
		'''
		_n_a 	 = self.n_axs(l, xc)
		_n_w 	 = self.n_ws(l)
		_La 	 = (_n_a**2 - 1) / (_n_a**2 + 2)
		_Lw 	 = (_n_w**2 - 1) / (_n_w**2 + 2)
		
		# Calculate the density of the dry air components, using the partial pressure of the dry air, i.e. p_dry = p_atm - p_water
		_rho_a 	 = self.rho(p=p-RH*self.svp(T), T=T, RH=0, xc=xc)				# Dry air component for the actual conditions.
		_rho_axs = self.rho(p=101325, T=288.15, RH=0, xc=xc) 					# Dry air component at standard conditions.
		
		if RH == 0:
			# nan is returned for _rho_w if we put in a (partial) pressure of 0
			_rho_w = 0
		else:
			# Partial pressure of water vapour is p_water = h*svp
			_rho_w 	 = self.rho(p=RH*self.svp(T), T=T, RH=RH, xc=xc)			# Moist air component for the actual conditions.
		_rho_ws	 = self.rho(p=1333, T=293.15, RH=1333/self.svp(293.15), xc=xc) 	# Moist air component at standard conditions.
		
		_Ltotal  = ( _rho_a / _rho_axs ) * _La + ( _rho_w / _rho_ws ) * _Lw

		# print 'Dry air component = ', (_rho_a/_rho_axs)*(_n_a-1)*1e8
		# print 'Wet component = ', ( _rho_w / _rho_ws ) * (_n_w - 1)*1e8
		# print 'Sum of components = ', (_rho_a/_rho_axs)*(_n_a-1)*1e8+( _rho_w / _rho_ws ) * (_n_w - 1)*1e8
		return np.sqrt( (1 + 2 * _Ltotal) / (1 - _Ltotal) )


	def dn_tph(self, l, T, p, RH, xc, dl=None, dT=None, dP=None, dRH=None, dCO2=None):
		'''
			Calculates the uncertainty in the refractive index, following the Ciddor model,
			by using fully propagated error analysis.
			This involves almost endless use of the chain rule and evaluation of both the moist
			air and dry air components. Note: dry air does not contain any water, i.e. H=0, x_w=0, etc.

			Parameters
			----------
			l 	: float
				Wavelength of light in microns
			p 	: float
				Atmospheric pressure in Pascal
			T 	: float
				Atmospheric temperature in Kelvin
			RH 	: float
				Atmospheric relative humidity (0<=RH<=1)
			xc 	: float
				CO2 density in parts per million
			dl 	: float 
				Uncertainty in the wavelength in microns
			dT 	: float (default: 0.2)
				Uncertainty in the atmospheric temperature in Kelvin
			dP 	: float (default: 20)
				Uncertainty in the atmospheric pressure in Pascal
			dRH : float (default: 0.02)
				Uncertainty in the atmospheric relative humidity
			dCO2: float (default: 20)
				Uncertainty in the atmospheric CO2 density in parts per million

			Returns
			-------
			dn 	: float
				Uncertainty in the refractive index of air
		'''

		if dP == None:
			dP = self.dP
		if dT == None:
			dT = self.dT
		if dRH == None:
			dRH = self.dH
		if dCO2 == None:
			dCO2 = self.dCO2

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
		# _rho_obs = self.rho(p, T, RH, xc)
		# _ref	 = refraction(lat, h, p, _rho_obs)
		# _alpha   = _ref.rc * np.sin(np.deg2rad(z)) / (_ref.rc + _ref.H)
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
		_n_a1 = self.n_axs(l, xc)
		_n_w1 = self.n_ws(l)

		_L1_a = (_n_a1**2 - 1) / (_n_a1**2 + 2)  # Dry air component of L
		_L1_w = (_n_w1**2 - 1) / (_n_w1**2 + 2)
		_L1   = ( _rho_a / _rho_axs ) * _L1_a + ( _rho_w / _rho_ws ) * _L1_w # Overall L
		_n1   = np.sqrt( (1 + 2 * _L1) / (1 - _L1) )

		# Partial derivatives for derivation to l1
		_dn_wdl  = _cf * 1e-8 * (-1 * _w1 * 2/l**3 - _w2 * 4/l**5 - _w3 * 6/l**7)
		_dLdn_w1 = _rho_w / _rho_ws * 6 * _n_w1 / (_n_w1**2 + 2)**2
		_dn_adl  = 1e-8 * (1 + 0.534e-6 * (xc - 450)) * ( - (2 * _k1 * l ) / ((_k0 * l**2 - 1)**2) - (2 * _k3 * l) / ((_k2 * l**2 - 1)**2) )
		_dLdn_a1 = _rho_a / _rho_axs * 6 * _n_a1 / (_n_a1**2 + 2)**2
		_dndL1   = ( ( (1 + 2*_L1) / (1 - _L1)**2 ) + (2 / (1-_L1)) ) / (2 * np.sqrt( (1 + 2*_L1) / (1 - _L1) ) ) 
		_dndl   = _dndL1 * ( _dLdn_a1 * _dn_adl + _dLdn_w1 * _dn_wdl)

		# # Partial derivatives for derivation to l2


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
		_dxw_wdf	 = 1
		_dfdT 		 = 2*_cc*_t
		_dxw_wdp_sv	 = self.f(_p_w, T) / _p_sv
		_drho_wdxw_w = _p_w * _Ma / (_Z_w * _R * T) * (_Mw / _Ma - 1)

		_tA =  _drho_adp_a * _dp_adp_sv * _dp_svdT +  _drho_adTZ * (_dTZdT_a + _dTZdxw_a * (_dxw_adf * _dfdT + _dxw_adp_sv * _dp_svdT) ) + _drho_adxw_a * ( _dxw_adf * _dfdT + _dxw_adp_sv * _dp_svdT )
		_tB =  _drho_wdp_w * _dp_wdp_sv * _dp_svdT +  _drho_wdTZ * (_dTZdT_w + _dTZdxw_w * (_dxw_wdf * _dfdT + _dxw_wdp_sv * _dp_svdT) ) + _drho_wdxw_w * ( _dxw_wdf * _dfdT + _dxw_wdp_sv * _dp_svdT )
		_dL1drho_a 	 = _L1_a / _rho_axs
		_dL1drho_w 	 = _L1_w / _rho_ws
		
		_dndT = _tA * (_dndL1 * _dL1drho_a) + _tB * (_dndL1 * _dL1drho_w)

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
		_dndP = _pA * (_dndL1 * _dL1drho_a)

		
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
		_dndRH = _hA * (_dndL1 * _dL1drho_a) + _hB * (_dndL1 * _dL1drho_w)


		# Partial derivatives for derivation to xc
		# Dry air parts:
		_drho_adMa   = _p_a * (1 - _x_w_a) / (_Z_a * _R * T)
		_dL1drho_axs = - (_rho_a * _L1_a)/(_rho_axs**2)
		_drho_axsdMa = _drho_adMa
		_dL1dL_a	 = _rho_a / _rho_axs
		_dL_adn_a1	 = 6 * _n_a1 / (_n_a1**2 + 2)**2

		# Moist air parts:
		_drho_wdMa   = _p_w * (1 - _x_w_w) / (_Z_w * _R * T)
		_dL1drho_ws = - (_rho_w * _L1_w)/(_rho_ws**2)
		_drho_wsdMa = _drho_wdMa

		# Other parts:
		_dMadx_c   = 12.011e-9
		_dn_adx_c1 = (self.n_as(l) - 1) * 0.534e-6

		_dL1dx_c = _dL1drho_a * _drho_adMa * _dMadx_c + _dL1drho_w * _drho_wdMa * _dMadx_c + _dL1drho_axs * _drho_axsdMa * _dMadx_c + _dL1drho_ws * _drho_wsdMa * _dMadx_c + _dL1dL_a * _dL_adn_a1 * _dn_adx_c1
		_dndCO2 = _dndL1 * _dL1dx_c

		# Final error:
		_nError = np.sqrt( (_dndl * dl)**2 + (_dndT * dT)**2  + (_dndP * dP)**2 + (_dndRH * dRH)**2 + (_dndCO2 * dCO2)**2 )

		return _nError	