#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Parabole Mirror Shift"""

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

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Thickness = 1000.0
P_Obj.Diameter = 300
P_Obj.Drawing = 0

# ______________________________________#

M1 = Kos.surf()
M1.Rc = -2 * P_Obj.Thickness
M1.Thickness = M1.Rc / 2
M1.k = -1.0
M1.Glass = "MIRROR"
M1.Diameter = 300
M1.ShiftY = 200

aa = 100
bb = 100


radio = 150

px = [radio * np.cos(np.radians(0)),
     radio * np.cos(np.radians(72)),
     radio * np.cos(np.radians(144)),
     radio * np.cos(np.radians(216)),
     radio * np.cos(np.radians(288)),
     radio * np.cos(np.radians(0))]

py = [radio * np.sin(np.radians(0)),
     radio * np.sin(np.radians(72)),
     radio * np.sin(np.radians(144)),
     radio * np.sin(np.radians(216)),
     radio * np.sin(np.radians(288)),
     radio * np.sin(np.radians(0))]

radio = 50

px1 = [radio * np.cos(np.radians(0)),
     radio * np.cos(np.radians(72)),
     radio * np.cos(np.radians(144)),
     radio * np.cos(np.radians(216)),
     radio * np.cos(np.radians(288)),
     radio * np.cos(np.radians(0))]

py1 = [radio * np.sin(np.radians(0)),
     radio * np.sin(np.radians(72)),
     radio * np.sin(np.radians(144)),
     radio * np.sin(np.radians(216)),
     radio * np.sin(np.radians(288)),
     radio * np.sin(np.radians(0))]

px.extend(px1)
py.extend(py1)



M1.UDA = [px, py]

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Glass = "AIR"
P_Ima.Diameter = 1600.0
P_Ima.Drawing = 0
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, M1, P_Ima]
configuracion_1 = Kos.Setup()

# ______________________________________#

Espejo = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Espejo)

# ______________________________________#



diametro = 300
num_puntos_lado = 10

x = np.linspace(-diametro/2, diametro/2, num_puntos_lado)
y = np.linspace(-diametro/2, diametro/2, num_puntos_lado)

X, Y = np.meshgrid(x, y)
X = X.ravel()
Y = Y.ravel()

for i in range(0,len(X)):
    x_0 = X[i]
    y_0 = Y[i]
    r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
    if r < diametro / 2.0:
        tet = 0.0
        pSource_0 = [x_0, y_0, 0.0]
        dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
        W = 0.4
        Espejo.Trace(pSource_0, dCos, W)
        Rayos.push()

# ______________________________________#

Kos.display3d(Espejo, Rayos, 0)



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

x,y,z,l,m,n = Rayos.pick(-1, coordinates="local")

print(R_RMS_delta(z, l, m, n, x, y))
