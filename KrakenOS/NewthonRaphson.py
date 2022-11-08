#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 18:07:46 2021

@author: joelherreravazquez
"""
import numpy as np
from scipy.optimize import fsolve

def f(z1):
    r=574
    L=0
    M=0
    N=1
    x0=10.0
    y0=0

    x1 = -z1*(L/N) - x0
    y1 = -z1*(M/N) - y0

    s= np.sqrt(x1**2 +y1**2)
    c= 1/r
    ze = (c*s**2)/np.sqrt(1-((c**2)*(s**2))) - r
    zt = np.abs(ze - z1)
    return(zt)

V=0
root = fsolve(f, V)
print(" Valor con fsolve", root)

z = 0
h= 0.0000001
for i in range(0,5):
    z = z - f(z)/( (f(z+h) - f(z-h) )/ (h*2.))
print("Valor con Newthon - Raphson: ",z)

