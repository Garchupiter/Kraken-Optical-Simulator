#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Perfect Lens"""

import time
import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("No instalado")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos

import scipy

def R_RMS_delta(Z1, L, M, N, X0, Y0):
    X1 = ((L / N) * Z1) + X0
    Y1 = ((M / N) * Z1) + Y0
    cenX = np.mean(X1)
    cenY = np.mean(Y1)
    x1 = (X1 - cenX)
    y1 = (Y1 - cenY)
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    return R_RMS

def BestFocus(X, Y, Z, L, M, N, system, mod=1):
    delta_Z = 0
    ZZ = (L, M, N, X, Y)
    v = scipy.optimize.fsolve(R_RMS_delta, delta_Z, args=ZZ)
    if mod ==1:
        system.SDT[-2].Thickness = system.SDT[-2].Thickness + v[0]
        system.SetData()
        system.SetSolid()
    return system, v[0]



P_Obj = Kos.surf(Thickness = 60, Diameter = 30.0)
L1 = Kos.surf(Thin_Lens = 25.0, Thickness =10.0, Diameter = 30.0)
L2 = Kos.surf(Thin_Lens = 15.0, Thickness = 50.0, Diameter = 30.0)
P_Ima = Kos.surf(Diameter = 40.0, Name = "i")

A = [P_Obj, L1, L2, P_Ima]
config_1 = Kos.Setup()

Lens = Kos.system(A, config_1)
Rays = Kos.raykeeper(Lens)

Surf, W, AperVal, AperType = 1, 0.45, L1.Diameter, "EPD"
P = Kos.PupilCalc(Lens, Surf, W, AperType, AperVal)
P.Samp, P.Ptype, P.FieldX, P.FieldType = 4, "fanx", 15, "height"

x, y, z, L, M, N = P.Pattern2Field()

Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)
x1,y1,z1,l1,m1,n1 = Rays.pick(-1, coordinates="local")

Lens, bf = BestFocus(x1, y1, z1, l1, m1, n1, Lens, mod=1)

Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)

x,y,z,l,m,n = Rays.pick(-1, coordinates="local")
print(x[0], y[0], z[0])
x,y,z,l,m,n = Rays.pick(-1, coordinates="global")
print(x[0], y[0], z[0])

Prx = Lens.Parax(W)



M, S_M, N_M, a, b, c, d, EFFL, PPA, PPP, C, N, D = Prx
print(EFFL, PPA, PPP)
# print(S_M)

Kos.display2d(Lens, Rays, 1, arrow=1)








