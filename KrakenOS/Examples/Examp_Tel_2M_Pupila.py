# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp  TEl 2M Pupila"""

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(1, '/Users/joelherreravazquez/Documents/GitHub/Kraken-Optical-Simulator/')

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
P_Obj.Rc = 0
P_Obj.Thickness = 1000 + 3.452200000000000E+003
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1.059E+003 * 2.0

# ______________________________________#

Thickness = 3.452200000000000E+003
M1 = Kos.surf()
M1.Rc = -9.638000000004009E+003
M1.Thickness = -Thickness
M1.k = -1.077310000000000E+000
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
M1.InDiameter = 250 * 2.0

# ______________________________________#

M2 = Kos.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness + 1.037525880125084E+003
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0
M2.TiltY = 0.1
M2.TiltX = 0.1
M2.AxisMove = 0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 300.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()
Telescopio = Kos.system(A, configuracion_1)

# ______________________________________#

W = 0.4
sup = 1
AperVal = 2010
AperType = "EPD"  # "STOP"
Pup = Kos.PupilCalc(Telescopio, sup, W, AperType, AperVal)

# ______________________________________#

print("Radio pupila de entrada: ")
print(Pup.RadPupInp)
print("Posicion pupila de entrada: ")
print(Pup.PosPupInp)
print("Radio pupila de salida: ")
print(Pup.RadPupOut)
print("Posicion pupila de salida: ")
print(Pup.PosPupOut)
print("Posicion pupila de salida respecto al plano focal: ")
print(Pup.PosPupOutFoc)
print("Orientación pupila de salida")
print(Pup.DirPupSal)
[L, M, N] = Pup.DirPupSal
print(L, M, N)
TetX = np.rad2deg(np.arcsin(-M))
TetY = np.rad2deg(np.arcsin(L / np.cos(np.arcsin(-M))))
print(TetX, TetY)
print("---------------------------------------------------------------")

# ______________________________________#

Pup.Samp = 10
Pup.Ptype = "hexapolar"
Pup.FieldY = 0.0
Pup.FieldType = "angle"
x, y, z, L, M, N = Pup.Pattern2Field()
Rayos = Kos.raykeeper(Telescopio)

# ______________________________________#

for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    W = 0.4
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()

# ______________________________________#

Kos.display2d(Telescopio, Rayos, 1, 1)

# ______________________________________#

X, Y, Z, L, M, N = Rayos.pick(-1)
plt.figure(300)
plt.plot(X, Y, 'x')
plt.axis('square')
plt.show(block=False)
