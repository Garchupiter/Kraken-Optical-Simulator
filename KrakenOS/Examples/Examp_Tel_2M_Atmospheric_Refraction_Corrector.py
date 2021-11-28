# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp-Tel_2M_Atmospheric_Refraction_Corrector"""

import pkg_resources
required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("No instalado")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos
import numpy as np
import matplotlib.pyplot as plt
import scipy

# _________________________________________________________________#
P_Obj = Kos.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 1000 + 3452.2
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1059.0 * 2.0
P_Obj.Drawing = 0
# _________________________________________________________________#
Thickness = 3452.2
M1 = Kos.surf()
M1.Rc = -9638.0
M1.Thickness = -Thickness
M1.k = -1.07731
M1.Glass = "MIRROR"
M1.Diameter = 1059.0 * 2.0
M1.InDiameter = 250 * 2.0
# _________________________________________________________________#
M2 = Kos.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness + 1037.525 - 300.0
M2.k = -4.3281
M2.Glass = "MIRROR"
M2.Diameter = 336.5 * 2.0
M2.AxisMove = 0
# _________________________________________________________________#
C1 = Kos.surf()
C1.Thickness = 5
C1.Glass = "BK7"
C1.Diameter = 100

C2 = Kos.surf()
C2.Thickness = C1.Thickness
C2.Glass = "F2"
C2.Diameter = 100
C2.TiltY = 1.55
C2.AxisMove = 0

C3 = Kos.surf()
C3.Thickness = C1.Thickness + 288.631
C3.Glass = "AIR"
C3.Diameter = 100
C3.TiltY = 0
C3.AxisMove = 0

# -------------------------------------

C6 = Kos.surf()
C6.Thickness = C1.Thickness
C6.Glass = "AIR"
C6.Diameter = 100

# _________________________________________________________________#
P_Ima = Kos.surf()
P_Ima.Diameter = 100.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"
A = [P_Obj, M1, M2, C1, C2, C3, P_Ima]
# _________________________________________________________________#
configuracion_1 = Kos.Setup()
Telescopio = Kos.system(A, configuracion_1)

Rayos1 = Kos.raykeeper(Telescopio)
Rayos2 = Kos.raykeeper(Telescopio)
Rayos3 = Kos.raykeeper(Telescopio)
Rayos = Kos.raykeeper(Telescopio)

W = 0.60169
sup = 1  # Difining M1 as enter pupil diameter
AperVal = 2000
AperType = "EPD"  # "STOP"
Pup = Kos.PupilCalc(Telescopio, sup, W, AperType, AperVal)
Pup.Samp = 5
Pup.FieldType = "angle"

Pup.AtmosRef = 1
Pup.T = 283.15  # k
Pup.P = 101300  # Pa
Pup.H = 0.5  # Humidity ratio 1 to 0
Pup.xc = 400  # ppm
Pup.lat = 31  # degrees
Pup.h = 2800  # meters
Pup.l1 = W  # micron
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


for i in range(0, len(xc)):
    pSource_0 = [xc[i], yc[i], zc[i]]
    dCos = [Lc[i], Mc[i], Nc[i]]
    Telescopio.Trace(pSource_0, dCos, W3)
    Rayos3.push()
    Rayos.push()

for i in range(0, len(xb)):
    pSource_0 = [xb[i], yb[i], zb[i]]
    dCos = [Lb[i], Mb[i], Nb[i]]
    Telescopio.Trace(pSource_0, dCos, W2)
    Rayos2.push()
    Rayos.push()

for i in range(0, len(xa)):
    pSource_0 = [xa[i], ya[i], za[i]]
    dCos = [La[i], Ma[i], Na[i]]
    Telescopio.Trace(pSource_0, dCos, W1)
    Rayos1.push()
    Rayos.push()

############################################

Kos.display3d(Telescopio, Rayos, 1)

X, Y, Z, L, M, N = Rayos1.pick(-1)
plt.plot(X * 1000.0, Y * 1000.0, 'x', c="b", ms=2)

X, Y, Z, L, M, N = Rayos3.pick(-1)
plt.plot(X * 1000.0, Y * 1000.0, 'x', c="g", ms=2)

X, Y, Z, L, M, N = Rayos2.pick(-1)
plt.plot(X * 1000.0, Y * 1000.0, 'x', c="r", ms=2)

# axis labeling
plt.xlabel('x')
plt.ylabel('y')

# figure name
plt.title('Spot Diagram')
plt.axis('square')
plt.ylim(-5, 5)
plt.show()

# v = Kos.BestFocus(X,Y,Z,L,M,N)
# rms = Kos.RMS(X,Y,Z,L,M,N)

# print("sol ---------------------")
# print(v)
