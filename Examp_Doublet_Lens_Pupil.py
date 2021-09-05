#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens Pupil"""

import KrakenOS as Kn

# _________________________________________#

P_Obj = Kn.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 100
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0
P_Obj.Name = "P_Obj"

# _________________________________________#

L1a = Kn.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# _________________________________________#

L1b = Kn.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# _________________________________________#

L1c = Kn.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001 - 40
L1c.Glass = "AIR"
L1c.Diameter = 30

# _________________________________________#

pupila = Kn.surf()
pupila.Rc = 30
pupila.Thickness = 40.
pupila.Glass = "AIR"
pupila.Diameter = 5
pupila.Name = "Pupila"
pupila.DespY = 0.
pupila.Nm_Poss=[-10,10]

# _________________________________________#

P_Ima = Kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 20.0
P_Ima.Name = "P_Ima"
P_Ima.Nm_Poss=[-10,10]

# _________________________________________#

A = [P_Obj, L1a, L1b, L1c, pupila, P_Ima]
config_1 = Kn.Kraken_setup()

# _________________________________________#

Doblete = Kn.system(A, config_1)
Rayos = Kn.raykeeper(Doblete)

# _________________________________________#

W = 0.4
sup = 4
AperVal = 10
AperType = "EPD"
Pup = Kn.PupilCalc(Doblete, sup, W, AperType, AperVal)

# _________________________________________#

print("Radio pupila de entrada: ")
print(Pup.RadPupInp)
print("Posición pupila de entrada: ")
print(Pup.PosPupInp)
print("Rádio pupila de salida: ")
print(Pup.RadPupOut)
print("Posicion pupila de salida: ")
print(Pup.PosPupOut)
print("Posicion pupila de salida respecto al plano focal: ")
print(Pup.PosPupOutFoc)
print("Orientación pupila de salida")
print(Pup.DirPupSal)

# _________________________________________#

[L, M, N] = Pup.DirPupSal
print(L, M, N)

# _________________________________________#

Pup.Samp = 3
Pup.Ptype = "fan"
Pup.FieldType = "angle"
Pup.FieldY = 2.0
x, y, z, L, M, N = Pup.Pattern2Field()
for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()

# _________________________________________#

Pup.FieldY = -2.0
x, y, z, L, M, N = Pup.Pattern2Field()
for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()

# _________________________________________#

Kn.display2d(Doblete, Rayos,0, 1)