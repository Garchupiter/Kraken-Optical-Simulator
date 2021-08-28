# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

import Kraken as Kn
import numpy as np
import matplotlib.pyplot as plt
import time

# _________________________________________________________________#
P_Obj = Kn.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 1000 + 3.452200000000000E+003
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1.059E+003 * 2.0
# _________________________________________________________________#
Thickness = 3.452200000000000E+003
M1 = Kn.surf()
M1.Rc = -9.638000000004009E+003
M1.Thickness = -Thickness
M1.k = -1.077310000000000E+000
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
M1.InDiameter = 250 * 2.0
# _________________________________________________________________#
M2 = Kn.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness + 1.037525880125084E+003
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0
M2.AxisMove = 0
# _________________________________________________________________#
P_Ima = Kn.surf()
P_Ima.Diameter = 1000.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"
A = [P_Obj, M1, M2, P_Ima]
# _________________________________________________________________#
configuracion_1 = Kn.Kraken_setup()
Telescopio = Kn.system(A, configuracion_1)
Rayos1 = Kn.raykeeper(Telescopio)
Rayos2 = Kn.raykeeper(Telescopio)
Rayos3 = Kn.raykeeper(Telescopio)
# _________________________________________________________________#
W = 0.4
sup = 1
AperVal = 2000
AperType = "EPD"  # "STOP"
Pup = Kn.PupilCalc(Telescopio, sup, W, AperType, AperVal)
Pup.Samp = 11
Pup.FieldType = "angle"

Pup.AtmosRef = 1
Pup.T = 283.15  # k
Pup.P = 101300  # Pa
Pup.H = 0.5  # Humidity ratio 1 to 0
Pup.xc = 400  # ppm
Pup.lat = 31  # degrees
Pup.h = 2800  # meters
Pup.l1 = 0.60169  # micron
Pup.l2 = 0.50169  # micron
Pup.z0 = 55.0  # degrees

Pup.Ptype = "hexapolar"
Pup.FieldX = 0.0

W1 = 0.50169
Pup.l2 = W1
xa, ya, za, La, Ma, Na = Pup.Pattern2Field()

W2 = 0.60169
Pup.l2 = W2
xb, yb, zb, Lb, Mb, Nb = Pup.Pattern2Field()

W3 = 0.70169
Pup.l2 = W3
xc, yc, zc, Lc, Mc, Nc = Pup.Pattern2Field()

############################################


for i in range(0, len(xa)):
    pSource_0 = [xa[i], ya[i], za[i]]
    dCos = [La[i], Ma[i], Na[i]]
    Telescopio.Trace(pSource_0, dCos, W1)
    Rayos1.push()

for i in range(0, len(xb)):
    pSource_0 = [xb[i], yb[i], zb[i]]
    dCos = [Lb[i], Mb[i], Nb[i]]
    Telescopio.Trace(pSource_0, dCos, W2)
    Rayos2.push()

for i in range(0, len(xc)):
    pSource_0 = [xc[i], yc[i], zc[i]]
    dCos = [Lc[i], Mc[i], Nc[i]]
    Telescopio.Trace(pSource_0, dCos, W3)
    Rayos3.push()

############################################

# Kn.display3d(Telescopio,Rayos,2)


X, Y, Z, L, M, N = Rayos1.pick(-1)
plt.plot(X * 1000.0, Y * 1000.0, 'x', c="b")

X, Y, Z, L, M, N = Rayos2.pick(-1)
plt.plot(X * 1000.0, Y * 1000.0, 'x', c="r")

X, Y, Z, L, M, N = Rayos3.pick(-1)
plt.plot(X * 1000.0, Y * 1000.0, 'x', c="g")

# axis labeling
plt.xlabel('x')
plt.ylabel('y')

# figure name
plt.title('Spot Diagram')
plt.axis('square')
plt.ylim(-np.pi, np.pi)
plt.show()
