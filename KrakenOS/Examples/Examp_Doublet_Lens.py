#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens"""


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


# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

L1a = Kos.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# ______________________________________#

L1b = Kos.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# ______________________________________#

L1c = Kos.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 3.0
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, L1a, L1b, L1c, P_Ima]
config_1 = Kos.Setup()

# ______________________________________#

Doblete = Kos.system(A, config_1)
Rayos1 = Kos.raykeeper(Doblete)
Rayos2 = Kos.raykeeper(Doblete)
Rayos3 = Kos.raykeeper(Doblete)
RayosT = Kos.raykeeper(Doblete)

# ______________________________________#

tam = 10
rad = 10.0
tsis = len(A) - 1
for j in range(-tam, tam + 1):
    for i in range(-tam, tam + 1):
        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            W = 0.4
            Doblete.Trace(pSource_0, dCos, W)
            Rayos1.push()
            RayosT.push()
            W = 0.5
            Doblete.Trace(pSource_0, dCos, W)
            Rayos2.push()
            RayosT.push()
            W = 0.6
            Doblete.Trace(pSource_0, dCos, W)
            Rayos3.push()
            RayosT.push()

# ______________________________________#

Kos.display2d(Doblete, RayosT, 0)

# ______________________________________#

X, Y, Z, L, M, N = RayosT.pick(-1)
dz = 0.0
for i in range(0, 10):
    FdeZ = R_RMS(L, M, N, X, Y, dz)
    derFdeZ = DER_R_RMS(L, M, N, X, Y, dz)
    print(dz, FdeZ, derFdeZ)
    dz = dz - (FdeZ / derFdeZ)
