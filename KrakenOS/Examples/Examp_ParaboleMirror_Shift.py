#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Parabole Mirror Shift"""

import numpy as np
import pkg_resources
import pickle


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
# M1.DerPres= 0.04

aa = 100
bb = 100

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

# Espejo = Kos.system_Lite(A, configuracion_1)





# with open('mi_objeto.pkl', 'wb') as archivo_salida:
#     # Usa pickle.dump para serializar y guardar el objeto en el archivo.
#     pickle.dump(Espejo, archivo_salida)

# with open('mi_objeto.pkl', 'rb') as archivo_entrada:
#     Espejo = pickle.load(archivo_entrada)



Rayos = Kos.raykeeper(Espejo)

# ______________________________________#

tam = 15
rad = 150.0
tsis = len(A) - 1
for i in range(-tam, tam + 1):
    for j in range(-tam, tam + 1):
        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            W = 0.4
            Espejo.Trace(pSource_0, dCos, W)
            Rayos.push()

# ______________________________________#

Kos.display2d(Espejo, Rayos, 0)



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
X, Y, Z, L, M, N = Rayos.pick(-1)

# ______________________________________#


import matplotlib.pyplot as plt
plt.plot(X, Y, 'x')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Spot Diagram')
plt.axis('square')
plt.show()
