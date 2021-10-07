#!/usr/bin/env Python3
import numpy as np

def sla_atms(RT, TT, DNT, GAMAL, R):
	'''
		Internal routine used by sla_refro(), ported from the SLALIB 
		FORTRAN code (Copyright P.T.Wallace.  All rights reserved.)

		The method is also described in "Explanatory Supplement 
		Astronomical Almanac" by P.K. Seidelmann, 1992, section 3.281.

		Refractive index and derivative with respect to height for the
		stratosphere.

		Parameters
		----------
		RT 		: 	float
			Height of tropopause from center of the Earth (meters)
		TT 		: 	float
			Temperature at the tropopause (K)
		DNT     : 	float
			Refractive index at the tropopause
		GAMAL 	: 	float
			Constant of the atmospheric model = G*MD/R
		R 		: 	float
			Current distance from the centre of the Earth (metre)

		Returns
		-------
		DN 		: 	float
			Refractive index at R
		RDNDR 	: 	float
			R * rate the refractive index is changing at R
	'''
	_B 		= GAMAL / TT
	_W 		= (DNT - 1) * np.exp( -_B * (R - RT) )
	DN  	= 1+_W
	RDNDR 	= -R * _B * _W
	return DN, RDNDR


def sla_atmt(R0, T0, ALPHA, GAMM2, DELM2, C1, C2, C3, C4, C5, C6, R):
	'''
		Internal routine used by sla_refro(), ported from the SLALIB 
		FORTRAN code (Copyright P.T.Wallace.  All rights reserved.)

		The method is also described in "Explanatory Supplement 
		Astronomical Almanac" by P.K. Seidelmann, 1992, section 3.281.

		Refractive index and derivative with respect to height for the
		troposphere.

		Parameters
		----------
		R0 		: 	float
			Height of observer from center of the Earth (meter)
		T0 		: 	float
			Temperature at the observer (K)
		ALPHA 	: 	float
			Alpha
		GAMM2 	: 	float 
			Gamma minus 2 - see HMNAO paper
		DELM2 	: 	float
			Delta minus 2
		C1  	: 	float
			Useful term
		C2  	: 	float
			Useful term
		C3  	: 	float
			Useful term
		C4  	: 	float
			Useful term
		C5  	: 	float
			Useful term
		C6  	: 	float
			Useful term
		R    	:  	float
			Current distance from the center of the Earth (meters)

		Returns
		-------
		T 		: 	float
			Temperature at R (K)
		DN 		:   float 
			Refractive index at R
		RDNDR 	: 	float
			R * rate the refractive index is changing at R

		Note that in the optical case C5 and C6 are zero.
	'''
	T 		= max(min(T0 - ALPHA * (R - R0), 320), 100)
	TT0 	= T / T0
	TT0GM2 	= TT0**GAMM2
	TT0DM2 	= TT0**DELM2
	DN 		= 1 + ( C1 * TT0GM2 - (C2 - C5 / T) * TT0DM2) * TT0
	RDNDR 	= R * (-C3 * TT0GM2 + (C4 - C6 / TT0) * TT0DM2)
	return T, DN, RDNDR


def sla_drange(a):
	'''
		Helper function for sla_refro().
	'''
	drange = a % (2*np.pi)
	if drange >= np.pi:
		drange -= np.abs(2*np.pi)*np.sign(a)
	return drange


def sla_refi(dn, rdndr):
	'''
		Helper function for sla_refro(), that calculates
		the refraction integrand.
	'''
	return rdndr/(dn+rdndr)

def sla_refro(ZOBS, HM, TDK, PMB, RH, WL, PHI, TLR, EPS):
	'''
		This is a direct port from the original FORTRAN SLA_REFRO()
		routine.

		Calculates Atmospheric refraction for radio and 
		optical/IR wavelengths, ported from the original SLALIB
		FORTRAN code (Copyright P.T.Wallace.  All rights reserved.).

		The method is also described in "Explanatory Supplement 
		Astronomical Almanac" by P.K. Seidelmann, 1992, section 3.281.
		
		Parameters
		----------
		ZOBS 	:	float
			Observed zenith distance of the source (radian)
		HM 		:	float
			Height of the observer above sea level (meter)
		TDK 	:	float
			Ambient temperature at the observer (K)
		PMB 	:	float
			Pressure at the observer (millibar)
		RH 		:	float
			Relative humidity at the observer (range 0-1)
		WL 		:	float
			Effective wavelength of the source (micrometre)
		PHI 	:	float
			Latitude of the observer (radian, astronomical)
		TLR 	:	float
			Temperature lapse rate in the troposphere (K/metre)
		EPS 	:	float
			Precision required to terminate iteration (radian)
		
		Returns
		-------
		REF 	:	float
			Refraction: in vacuo ZD minus observed ZD (radian)
		
		Notes
		-----
		1  A suggested value for the TLR argument is 0.0065D0.  The
		   refraction is significantly affected by TLR, and if studies
		   of the local atmosphere have been carried out a better TLR
		   value may be available.  The sign of the supplied TLR value
		   is ignored.
		2  A suggested value for the EPS argument is 1D-8.  The result is
		   usually at least two orders of magnitude more computationally
		   precise than the supplied EPS value.
		3  The routine computes the refraction for zenith distances up
		   to and a little beyond 90 deg using the method of Hohenkerk
		   and Sinclair (NAO Technical Notes 59 and 63, subsequently adopted
		   in the Explanatory Supplement, 1992 edition - see section 3.281).
		4  The code is a development of the optical/IR refraction subroutine
		   AREF of C.Hohenkerk (HMNAO, September 1984), with extensions to
		   support the radio case.  Apart from merely cosmetic changes, the
		   following modifications to the original HMNAO optical/IR refraction
		   code have been made:
		   .  The angle arguments have been changed to radians.
		   .  Any value of ZOBS is allowed (see note 6, below).
		   .  Other argument values have been limited to safe values.
		   .  Murray's values for the gas constants have been used
		      (Vectorial Astrometry, Adam Hilger, 1983).
		   .  The numerical integration phase has been rearranged for
		      extra clarity.
		   .  A better model for Ps(T) has been adopted (taken from
		      Gill, Atmosphere-Ocean Dynamics, Academic Press, 1982).
		   .  More accurate expressions for Pwo have been adopted
		      (again from Gill 1982).
		   .  The formula for the water vapour pressure, given the
		      saturation pressure and the relative humidity, is from
		      Crane (1976), expression 2.5.5.
		   .  Provision for radio wavelengths has been added using
		      expressions devised by A.T.Sinclair, RGO (private
		      communication 1989).  The refractivity model currently
		      used is from J.M.Rueger, "Refractive Index Formulae for
		      Electronic Distance Measurement with Radio and Millimetre
		      Waves", in Unisurv Report S-68 (2002), School of Surveying
		      and Spatial Information Systems, University of New South
		      Wales, Sydney, Australia.
		   .  The optical refractivity for dry air is from Resolution 3 of
		      the International Association of Geodesy adopted at the XXIIth
		      General Assembly in Birmingham, UK, 1999.
		   .  Various small changes have been made to gain speed.
		5  The radio refraction is chosen by specifying WL > 100 micrometres.
		   Because the algorithm takes no account of the ionosphere, the
		   accuracy deteriorates at low frequencies, below about 30 MHz.
		6  Before use, the value of ZOBS is expressed in the range +/- pi.
		   If this ranged ZOBS is -ve, the result REF is computed from its
		   absolute value before being made -ve to match.  In addition, if
		   it has an absolute value greater than 93 deg, a fixed REF value
		   equal to the result for ZOBS = 93 deg is returned, appropriately
		   signed.
		7  As in the original Hohenkerk and Sinclair algorithm, fixed values
		   of the water vapour polytrope exponent, the height of the
		   tropopause, and the height at which refraction is negligible are
		   used.
		8  The radio refraction has been tested against work done by
		   Iain Coulson, JACH, (private communication 1995) for the
		   James Clerk Maxwell Telescope, Mauna Kea.  For typical conditions,
		   agreement at the 0.1 arcsec level is achieved for moderate ZD,
		   worsening to perhaps 0.5-1.0 arcsec at ZD 80 deg.  At hot and
		   humid sea-level sites the accuracy will not be as good.
		9  It should be noted that the relative humidity RH is formally
		   defined in terms of "mixing ratio" rather than pressures or
		   densities as is often stated.  It is the mass of water per unit
		   mass of dry air divided by that for saturated air at the same
		   temperature and pressure (see Gill 1982).
		10 The algorithm is designed for observers in the troposphere.  The
		   supplied temperature, pressure and lapse rate are assumed to be
		   for a point in the troposphere and are used to define a model
		   atmosphere with the tropopause at 11km altitude and a constant
		   temperature above that.  However, in practice, the refraction
		   values returned for stratospheric observers, at altitudes up to
		   25km, are quite usable.
	'''

	# Fixed parameters
	D93 	= np.deg2rad(93)# 93 degrees in radians
	GCR 	= 8314.32 		# Universal gas constant
	DMD 	= 28.9644 		# Molecular weight of dry air
	DMW 	= 18.0152 		# Molecular weight of water vapour
	S 		= 6378120 		# Mean Earth radius (metre)
	DELTA 	= 18.36 		# Exponent of temperature dependence of water vapour pressure
	HT 		= 11000			# Height of tropopause (metre)
	HS 		= 80000			# Upper limit for refractive effects (metre)
	ISMAX 	= 16384			# Numerical integration: maximum number of strips.
	
	# Transform ZOBS into the normal range.
	ZOBS1 = sla_drange(ZOBS)
	ZOBS2 = min(abs(ZOBS1),D93)

	# Keep other arguments within safe bounds.
	HMOK  = min(max(HM,-1E3),HS)
	TDKOK = min(max(TDK,100),500)
	PMBOK = min(max(PMB,0),10000)
	RHOK  = min(max(RH,0),1)
	WLOK  = max(WL,0.1)
	ALPHA = min(max(abs(TLR),0.001),0.01)

	# Tolerance for iteration.
	TOL   = min(max(abs(EPS),1E-12),0.1)/2

	# Decide whether optical/IR or radio case - switch at 100 microns.
	OPTIC = (WLOK<=100)

	# Set up model atmosphere parameters defined at the observer.
	WLSQ  = WLOK*WLOK

	# Local gravity, can be made consistent with AstroAtmosphere
	GB    = 9.784*(1-0.0026*np.cos(PHI+PHI)-0.00000028*HMOK)

	if OPTIC:
		# Refractivity from Barrell & Sears (1939)!
		A = (287.6155+(1.62887+0.01360/WLSQ)/WLSQ)*273.15E-6/1013.25
	else:
		A = 77.6890E-6

	GAMAL = (GB*DMD)/GCR
	GAMMA = GAMAL/ALPHA
	GAMM2 = GAMMA-2
	DELM2 = DELTA-2
	TDC   = TDKOK-273.15
	PSAT  = 10**((0.7859+0.03477*TDC)/(1+0.00412*TDC))\
				*(1+PMBOK*(4.5E-6+6E-10*TDC*TDC))

	if PMBOK > 0:
		PWO = RHOK*PSAT/(1-(1-RHOK)*PSAT/PMBOK)
	else:
		PWO = 0

	W  = PWO*(1-DMW/DMD)*GAMMA/(DELTA-GAMMA)
	C1 = A*(PMBOK+W)/TDKOK
	if OPTIC:
		C2 = (A*W+11.2684E-6*PWO)/TDKOK
	else:
		C2 = (A*W+6.3938E-6*PWO)/TDKOK
	
	C3 = (GAMMA-1)*ALPHA*C1/TDKOK
	C4 = (DELTA-1)*ALPHA*C2/TDKOK
	if OPTIC:
		C5 = 0
		C6 = 0
	else:
		C5 = 375463E-6*PWO/TDKOK
		C6 = C5*DELM2*ALPHA/(TDKOK*TDKOK)

	# Conditions at the observer.
	R0  = S+HMOK
	TEMPO, DN0, RDNDR0 = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,R0)
	SK0 = DN0*R0*np.sin(ZOBS2)
	F0  = sla_refi(DN0,RDNDR0)

	# Conditions in the troposphere at the tropopause.
	RT = S+max(HT,HMOK)
	TT, DNT, RDNDRT = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,RT)
	SINE = SK0/(RT*DNT)
	ZT = np.arctan2(SINE,np.sqrt(max(1-SINE*SINE,0)))
	FT = sla_refi(DNT,RDNDRT)

	# Conditions in the stratosphere at the tropopause.
	DNTS,RDNDRP = sla_atms(RT,TT,DNT,GAMAL,RT)
	SINE = SK0/(RT*DNTS)
	ZTS  = np.arctan2(SINE,np.sqrt(max(1-SINE*SINE,0)))
	FTS  = sla_refi(DNTS,RDNDRP)

	# Conditions at the stratosphere limit.
	RS   = S+HS
	DNS,RDNDRS = sla_atms(RT,TT,DNT,GAMAL,RS)
	SINE = SK0/(RS*DNS)
	ZS   = np.arctan2(SINE,np.sqrt(max(1-SINE*SINE,0)))
	FS   = sla_refi(DNS,RDNDRS)

	# Variable initialization to avoid compiler warning.
	REFT = 0

	# Integrate the refraction integral in two parts;  first in the
	# troposphere (K=1), then in the stratosphere (K=2).

	for K in [1,2]:
		# Initialize previous refraction to ensure at least two iterations.
		REFOLD = 1

		# Start off with 8 strips.
		IS = 8

		# Start Z, Z range, and start and end values.
		if K==1:
			Z0 = ZOBS2
			ZRANGE = ZT-Z0
			FB = F0
			FF = FT
		else:
			Z0 = ZTS
			ZRANGE = ZS-Z0
			FB = FTS
			FF = FS

		# Sums of odd and even values.
		FO = 0
		FE = 0

		# First time through the loop we have to do every point.
		N = 1

		# Start of iteration loop (terminates at specified precision).
		LOOP = True
		while LOOP:

			# Strip width.
			H = ZRANGE/IS

			# Initialize distance from Earth centre for quadrature pass.
			if K==1:
				R = R0
			else:
				R = RT

			# One pass (no need to compute evens after first time).
			# DO I=1,IS-1,N
			for I in np.arange(1, IS-1, N):

				# Sine of observed zenith distance.
				SZ = np.sin(Z0+H*I)

				# Find R (to the nearest metre, maximum four iterations).
				if (SZ>1E-20):
					W  = SK0/SZ
					RG = R
					DR = 1E6
					J  = 0
					while (abs(DR)>1 and J<4):
						J+=1
						if K==1:
							TG,DN,RDNDR = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,RG)
						else:
							DN,RDNDR = sla_atms(RT,TT,DNT,GAMAL,RG)
						DR = (RG*DN-W)/(DN+RDNDR)
						RG = RG-DR
					R = RG

				# Find the refractive index and integrand at R.
				if K==1:
					T,DN,RDNDR = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,R)
				else:
					DN,RDNDR = sla_atms(RT,TT,DNT,GAMAL,R)

				F = sla_refi(DN,RDNDR)

				# Accumulate odd and (first time only) even values.
				if N==1 and (I%2)==0:
					FE = FE+F
				else:
					FO = FO+F


			# Evaluate the integrand using Simpson's Rule.
			REFP = H*(FB+4*FO+2*FE+FF)/3

			# Has the required precision been achieved (or can't be)?
			if (abs(REFP-REFOLD)>TOL and IS<ISMAX):
				# No: prepare for next iteration.

				# Save current value for convergence test.
				REFOLD = REFP

				# Double the number of strips.
				IS = IS+IS

				# Sum of all current values = sum of next pass's even values.
				FE = FE+FO

				# Prepare for new odd values.
				FO = 0

				# Skip even values next time.
				N = 2
			else:
				# Yes: save troposphere component and terminate the loop.
				if (K==1):
					REFT = REFP
				LOOP = False

	# Result.
	REF = REFT+REFP
	if (ZOBS1<0): 
		REF = -REF
	return REF


# =============================================================================
# WORK IN PROGRESS; Implementation of the slalib model for 
# 					any refractive index model.
# Main issue is in the calculation of dn/dr for any model
# of the refractive index. One possibility is following 
# the Corbard error function model, which assumes the
# atmospheric density follows the equation rho(r)=exp(-u/a).
# But then the SLA_REFRO routine would not have to be modified
# and one would simply use the refraction.corbard() function.
# =============================================================================
#
# def get_rc(lat):
# 		'''
# 			Returns the radius of the curvature at the observer in km,
# 			taking into account the latitude of the observer and the corresponding
# 			ellipsoidality of the earth.

# 			Returns
# 			-------
# 			rc 	: float
# 				Local radius of curvature of the Earth in kilometers.
# 		'''
# 		_A   = 6378.137 # km
# 		_B   = 6356.752 # km
# 		_phi = np.deg2rad(lat)
		
# 		_rc0 = (_A*_B)**2 / (_A**2 * np.cos(_phi)**2 + _B**2 * np.sin(_phi)**2)**(3/2)
# 		return _rc0*1000


# def get_local_g(lat, h):
# 	'''
# 		Returns the local gravitic acceleration for a given latitude 
# 		and altitude.

# 		Parameters
# 		----------
# 		lat : float
# 			Latitude in degrees
# 		h 	: float
# 			Altitude above sea level in meters

# 		Returns
# 		-------
# 		g 	: float
# 			Local acceleration due to gravity in m/s^2
# 	'''
# 	_phi = np.deg2rad(lat)
# 	_c1 = 5.2790414e-3
# 	_c2 = 2.32718e-5
# 	_c3 = 1.262e-7
# 	_c4 = 7e-10
# 	_g0 = 9.780327 # m/s^2
# 	_g0_local = _g0 * (1 + _c1 * np.sin(_phi)**2 + _c2 * np.sin(_phi)**4 + _c3 * np.sin(_phi)**6 + _c4 * np.sin(_phi)**8)
# 	_g = _g0_local - (3.0877e-6 - 4.3e-9 * np.sin(_phi)**2) * h + 7.2e-13 * h**2
# 	return _g

# def aa_refro(zenith, l, h, T, p, RH, xc, lat, TLR, EPS):
# 	'''
# 		Modified version of sla_refro(), updating some constants
# 		and changing the units on input and output parameters.

# 		The function calculates the atmospheric refraction
# 		optical/IR wavelengths (not radio!), ported from the original 
# 		SLALIB FORTRAN code (Copyright P.T.Wallace.  All rights 
# 		reserved.).

# 		The method is also described in "Explanatory Supplement 
# 		Astronomical Almanac" by P.K. Seidelmann, 1992, section 3.281.

# 		See also the description of sla_refro().
		
# 		Parameters
# 		----------
# 		zenith 	:	float
# 			Observed zenith distance of the source (degrees)
# 		l 		:	float
# 			Effective wavelength of the source (micrometre)
# 		h 		:	float
# 			Height of the observer above sea level (meter)
# 		T 		:	float
# 			Ambient temperature at the observer (K)
# 		p 		:	float
# 			Pressure at the observer (Pascal)
# 		RH 		:	float
# 			Relative humidity at the observer (range 0-1)
# 		xc 		: float
# 			CO2 density in parts per million
# 		lat 	:	float
# 			Latitude of the observer (degrees)
# 		TLR 	:	float
# 			Temperature lapse rate in the troposphere (K/metre)
# 		EPS 	:	float
# 			Precision required to terminate iteration (radian)
		
# 		Returns
# 		-------
# 		REF 	:	float
# 			Refraction: in vacuo ZD minus observed ZD (radian)
# 	'''
	
# 	# Converting some units for use later
# 	ZOBS = np.deg2rad(zenith)
# 	HM 	 = h
# 	TDK  = T
# 	PMB  = p / 100
# 	WL   = l
# 	PHI  = np.deg2rad(lat)

# 	# Fixed parameters
# 	D93 	= np.deg2rad(93)# 93 degrees in radians
# 	GCR 	= 8314.510 		# Universal gas constant
# 	DMD 	= 28.9635 + 12.011e-6 * (xc - 400) # Molar mass in g/mol of dry air containing xc ppm of CO2
# 	DMW 	= 18.0152 		# Molecular weight of water vapour
# 	S 		= get_rc(lat)	# Mean Earth radius (meters)
# 	DELTA 	= 18.36 		# Exponent of temperature dependence of water vapour pressure
# 	HT 		= 11000			# Height of tropopause (metre)
# 	HS 		= 80000			# Upper limit for refractive effects (meters)
# 	ISMAX 	= 16384			# Numerical integration: maximum number of strips.
	
# 	# Transform ZOBS into the normal range.
# 	ZOBS1 = sla_drange(ZOBS)
# 	ZOBS2 = min(abs(ZOBS1),D93)

# 	# Keep other arguments within safe bounds.
# 	HMOK  = min(max(HM,-1E3),HS)
# 	TDKOK = min(max(TDK,100),500)
# 	PMBOK = min(max(PMB,0),10000)
# 	RHOK  = min(max(RH,0),1)
# 	WLOK  = max(WL,0.1)
# 	ALPHA = min(max(abs(TLR),0.001),0.01)

# 	# Tolerance for iteration.
# 	TOL   = min(max(abs(EPS),1E-12),0.1)/2

# 	# Set up model atmosphere parameters defined at the observer.
# 	WLSQ  = WLOK*WLOK

# 	# Local gravity
# 	GB = get_local_g(lat, h)

# 	# Refractivity of air
# 	A = (287.6155+(1.62887+0.01360/WLSQ)/WLSQ)*273.15E-6/1013.25

# 	GAMAL = (GB * DMD) / GCR
# 	GAMMA = GAMAL / ALPHA
# 	GAMM2 = GAMMA - 2
# 	DELM2 = DELTA - 2
# 	TDC   = TDKOK - 273.15
# 	PSAT  = 10**( (0.7859+0.03477*TDC) / (1+0.00412*TDC) )\
# 				* (1 + PMBOK * (4.5E-6 + 6E-10 * TDC*TDC) )

# 	if PMBOK > 0:
# 		PWO = RHOK*PSAT/(1-(1-RHOK)*PSAT/PMBOK)
# 	else:
# 		PWO = 0

# 	W  = PWO*(1-DMW/DMD)*GAMMA/(DELTA-GAMMA)
# 	C1 = A*(PMBOK+W)/TDKOK
# 	C2 = (A*W+11.2684E-6*PWO)/TDKOK
	
# 	C3 = (GAMMA-1)*ALPHA*C1/TDKOK
# 	C4 = (DELTA-1)*ALPHA*C2/TDKOK
# 	C5 = 0
# 	C6 = 0

# 	# Conditions at the observer.
# 	R0  = S+HMOK
# 	TEMPO, DN0, RDNDR0 = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,R0)
# 	SK0 = DN0*R0*np.sin(ZOBS2)
# 	F0  = sla_refi(DN0,RDNDR0)

# 	# Conditions in the troposphere at the tropopause.
# 	RT = S+max(HT,HMOK)
# 	TT, DNT, RDNDRT = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,RT)
# 	SINE = SK0/(RT*DNT)
# 	ZT = np.arctan2(SINE,np.sqrt(max(1-SINE*SINE,0)))
# 	FT = sla_refi(DNT,RDNDRT)

# 	# Conditions in the stratosphere at the tropopause.
# 	DNTS,RDNDRP = sla_atms(RT,TT,DNT,GAMAL,RT)
# 	SINE = SK0/(RT*DNTS)
# 	ZTS  = np.arctan2(SINE,np.sqrt(max(1-SINE*SINE,0)))
# 	FTS  = sla_refi(DNTS,RDNDRP)

# 	# Conditions at the stratosphere limit.
# 	RS   = S+HS
# 	DNS,RDNDRS = sla_atms(RT,TT,DNT,GAMAL,RS)
# 	SINE = SK0/(RS*DNS)
# 	ZS   = np.arctan2(SINE,np.sqrt(max(1-SINE*SINE,0)))
# 	FS   = sla_refi(DNS,RDNDRS)

# 	# Variable initialization to avoid compiler warning.
# 	REFT = 0

# 	# Integrate the refraction integral in two parts;  first in the
# 	# troposphere (K=1), then in the stratosphere (K=2).
# 	for K in [1,2]:
# 		# Initialize previous refraction to ensure at least two iterations.
# 		REFOLD = 1

# 		# Start off with 8 strips.
# 		IS = 8

# 		# Start Z, Z range, and start and end values.
# 		if K==1:
# 			Z0 = ZOBS2
# 			ZRANGE = ZT-Z0
# 			FB = F0
# 			FF = FT
# 		else:
# 			Z0 = ZTS
# 			ZRANGE = ZS-Z0
# 			FB = FTS
# 			FF = FS

# 		# Sums of odd and even values.
# 		FO = 0
# 		FE = 0

# 		# First time through the loop we have to do every point.
# 		N = 1

# 		# Start of iteration loop (terminates at specified precision).
# 		LOOP = True
# 		while LOOP:

# 			# Strip width.
# 			H = ZRANGE/IS

# 			# Initialize distance from Earth centre for quadrature pass.
# 			if K==1:
# 				R = R0
# 			else:
# 				R = RT

# 			# One pass (no need to compute evens after first time).
# 			# DO I=1,IS-1,N
# 			for I in np.arange(1, IS-1, N):

# 				# Sine of observed zenith distance.
# 				SZ = np.sin(Z0+H*I)

# 				# Find R (to the nearest metre, maximum four iterations).
# 				if (SZ>1E-20):
# 					W  = SK0/SZ
# 					RG = R
# 					DR = 1E6
# 					J  = 0
# 					while (abs(DR)>1 and J<4):
# 						J+=1
# 						if K==1:
# 							TG,DN,RDNDR = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,RG)
# 						else:
# 							DN,RDNDR = sla_atms(RT,TT,DNT,GAMAL,RG)
# 						DR = (RG*DN-W)/(DN+RDNDR)
# 						RG = RG-DR
# 					R = RG

# 				# Find the refractive index and integrand at R.
# 				if K==1:
# 					T,DN,RDNDR = sla_atmt(R0,TDKOK,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,R)
# 				else:
# 					DN,RDNDR = sla_atms(RT,TT,DNT,GAMAL,R)

# 				F = sla_refi(DN,RDNDR)

# 				# Accumulate odd and (first time only) even values.
# 				if N==1 and (I%2)==0:
# 					FE = FE+F
# 				else:
# 					FO = FO+F


# 			# Evaluate the integrand using Simpson's Rule.
# 			REFP = H*(FB+4*FO+2*FE+FF)/3

# 			# Has the required precision been achieved (or can't be)?
# 			if (abs(REFP-REFOLD)>TOL and IS<ISMAX):
# 				# No: prepare for next iteration.

# 				# Save current value for convergence test.
# 				REFOLD = REFP

# 				# Double the number of strips.
# 				IS = IS+IS

# 				# Sum of all current values = sum of next pass's even values.
# 				FE = FE+FO

# 				# Prepare for new odd values.
# 				FO = 0

# 				# Skip even values next time.
# 				N = 2
# 			else:
# 				# Yes: save troposphere component and terminate the loop.
# 				if (K==1):
# 					REFT = REFP
# 				LOOP = False

# 	# Result.
# 	REF = REFT+REFP
# 	if (ZOBS1<0): 
# 		REF = -REF
# 	return REF



