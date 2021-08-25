#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 14:22:55 2021

@author: joelherreravazquez
"""

from AstroAtmosphere import *


Z0=75.0	    #	The observed zenith distance of the object in degrees
H0=0.0	    #	The height of the observer above sea level in metres
T0=283.15	#	The temperature at the observer in degrees K
P0=1005.0	#	The pressure at the observer in millibars
UPS=0.0	    #	The relative humidity at the observer
WL=0.50169	#	The wavelength of the light at the observer in micrometres
PH=50.0  	#	The latitude of the observer in degrees
AS=0.005694 #	The temperature lapse rate in degrees K/metre in the troposphere, the absolute value is used




# T - Temperature [K]
# p - Pressure [Pa]
# RH - Relative humidity
# xc - CO2 density [ppm]
# lat - Latitude [deg]
# h - Height above sea level [m]
# l - Wavelength(s) [um]
# z - Zenith angle [deg]


# Parameters at Cerro Armazones
T   = 283.15    # k
P   = 100500     # Pa
H   = 0.0
xc  = 450       # ppm
lat = 50    # degrees
h   = 0      # m




l1  = 5000.60169      # micron
l2  = 0.50169      # micron

for i in range(0,1000):
    L=l2+.1*i/1000.0
    z0  = 75.0
    
    # Initializing dispersion model
    at  = Observatory()
    
    # Calculating indices of refraction for l1 and l2
    n1  = at.n_tph(l=l1, T=T, p=P, RH=H, xc=xc)
    n2  = at.n_tph(l=L, T=T, p=P, RH=H, xc=xc)
    
    
    # Density of the atmosphere (following CIPM-81/91 equations)
    rho = at.rho(p=P, T=T, RH=H, xc=xc)
    
    # Initializing refraction model and setting the reduced height
    disp = dispersion(lat, h)
    disp.setReducedHeight(P, rho)
    
    # Calculation of the atmopheric dipsersion
    atm_dispersion = disp.cassini(n1, n2, z0) * 3.6e6
    print ('The dispersion is %.03f milli arc seconds' %(atm_dispersion), L)



