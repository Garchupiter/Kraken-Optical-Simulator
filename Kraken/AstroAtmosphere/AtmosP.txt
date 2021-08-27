#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 14:22:55 2021

@author: joelherreravazquez
"""

from AstroAtmosphere import *

# Parameters at Cerro Armazones
T   = 279.65    # k
P   = 71200     # Pa
H   = 0.22
xc  = 450       # ppm
lat = -24.5983  # degrees
h   = 3064      # m
l1  = 1.49      # micron
l2  = 1.78      # micron
z0  = 30

# Initializing dispersion model
at  = Observatory()

# Calculating indices of refraction for l1 and l2
n1  = at.n_tph(l=l1, T=T, p=P, RH=H, xc=xc)
n2  = at.n_tph(l=l2, T=T, p=P, RH=H, xc=xc)

# Density of the atmosphere (following CIPM-81/91 equations)
rho = at.rho(p=P, T=T, RH=H, xc=xc)

# Initializing refraction model and setting the reduced height
disp = dispersion(lat, h)
disp.setReducedHeight(P, rho)

# Calculation of the atmopheric dipsersion
atm_dispersion = disp.cassini(n1, n2, z0) * 3.6e6
print ('The dispersion is %.03f milli arc seconds' %(atm_dispersion))