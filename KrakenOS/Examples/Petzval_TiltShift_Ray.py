#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens"""

import time
import numpy as np

import pkg_resources
""" Looking for if KrakenOS is installed, if not, it assumes that
an folder downloaded from github is run"""

required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("Not installed")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos

# ______________________________________#

start_time = time.time()


# ______________________________________#

def R_RMS(L, M, N, X, Y, delta_Z):
    cenX = np.mean(X)
    cenY = np.mean(Y)
    nX = X - cenX
    nY = Y - cenY
    x1 = ((L / N) * delta_Z) + nX
    y1 = ((M / N) * delta_Z) + nY
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    return R_RMS


# ______________________________________#

def DER_R_RMS(L, M, N, X, Y, delta_Z):
    h = 0.001
    f1 = R_RMS(L, M, N, X, Y, delta_Z + h)
    f2 = R_RMS(L, M, N, X, Y, delta_Z - h)
    der = (f1 - f2) / (2.0 * h)
    return der






OBJ=Kos.surf()
OBJ.Rc=0.0
OBJ.Thickness=10.0
OBJ.Glass="AIR"
OBJ.Diameter=56


STO=Kos.surf()
STO.Rc=0.0
STO.Thickness=10.0
STO.Glass="AIR"
STO.Diameter=56


s2=Kos.surf()
s2.Rc=91.66860362426125
s2.Thickness=2.16760007E+1
s2.Glass="BALKN3"
s2.Diameter=56
s2.DespX = 1.0
s2.AxisMove = 0

s3=Kos.surf()
s3.Rc=-60.18807811401379
s3.Thickness=3.48399646
s3.Glass="F4"
s3.Diameter=56
s3.DespY = 2.0
s3.AxisMove = 0

s4=Kos.surf()
s4.Rc=3972.618614963046
s4.Thickness=7.65093542E+1
s4.Glass="AIR"
s4.Diameter=56
s4.TiltX = 3.0
s4.AxisMove = 0

s5=Kos.surf()
s5.Rc=35.206739696593594
s5.Thickness=1.7360596E+1
s5.Glass="BK7"
s5.Diameter=48
s5.TiltY = 4.0
s5.AxisMove = 0

s6=Kos.surf()
s6.Rc=-58.825151255787866
s6.Thickness=3.48399646
s6.Glass="F2"
s6.Diameter=48
s6.DespX = -1
s6.AxisMove = 0

s7=Kos.surf()
s7.Rc=-166.14088338204223
s7.Thickness=1.93995258E+1
s7.Glass="AIR"
s7.Diameter=48
s7.DespY = -2
s7.AxisMove = 0


s8=Kos.surf()
s8.Rc=-27.78333066229023
s8.Thickness=1.96964573
s8.Glass="F2"
s8.Diameter = 32
s8.TiltX = -3
s8.AxisMove = 0


s9=Kos.surf()
s9.Rc=1111723.0516534068
s9.Thickness=15.1962995812
s9.Glass="AIR"
s9.Diameter=32

IMA=Kos.surf()
IMA.Rc=0.0
IMA.Thickness=0.0
IMA.Glass="AIR"
IMA.Diameter=13.98677

A = [OBJ, STO, s2, s3, s4, s5, s6, s7, s8, s9, IMA]

config_1 = Kos.Setup()

Petzval = Kos.system(A, config_1)
Rays = Kos.raykeeper(Petzval)

# ______________________________________#

W=0.4861
pSource_0 = [0, 0, 0]
dCos = [0., 0., 1.]
Petzval.Trace(pSource_0, dCos, W)
print(1.0)
for i in range(1,11):
    P=Petzval.OST_XYZ[-i]
    print(P[0])

print(" ")
for i in range(1,11):
    P=Petzval.OST_XYZ[-i]
    print(P[1])


print(type(P[1]))
