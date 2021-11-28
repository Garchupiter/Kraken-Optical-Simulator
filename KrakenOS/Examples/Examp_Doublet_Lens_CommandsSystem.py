#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens Commands System"""

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

# _________________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 0.1
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.

# _________________________________________#

P_Obj2 = Kos.surf()
P_Obj2.Rc = 0.0
P_Obj2.Thickness = 10
P_Obj2.Glass = "AIR"
P_Obj2.Diameter = 30.0

# _________________________________________#

L1a = Kos.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# _________________________________________#

L1b = Kos.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# _________________________________________#

L1c = Kos.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

# _________________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 18.0
P_Ima.Name = "Plano imagen"

# _________________________________________#

A = [P_Obj, P_Obj2, L1a, L1b, L1c, P_Ima]
configuracion_1 = Kos.Setup()

# _________________________________________#

Doblete = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Doblete)

# _________________________________________#

pSource_0 = [0, 14, 0]
tet = 0.1
dCos = [0.0, np.sin(np.deg2rad(tet)), -np.cos(np.deg2rad(tet))]
W = 0.4
Doblete.Trace(pSource_0, dCos, W)
Rayos.push()

# _________________________________________#

Kos.display3d(Doblete, Rayos, 2)

# _________________________________________#

print("Distancia focal efectiva")
print(Doblete.EFFL)
print("Plano principal anterior")
print(Doblete.PPA)
print("Plano principal posterior")
print(Doblete.PPP)
print("Superficies tocadas por el rayo")
print(Doblete.SURFACE)
print("Nombre de la superficie")
print(Doblete.NAME)
print("Vidrio de la superficie")
print(Doblete.GLASS)
print("Coordenadas del rayo en las superficies")
print(Doblete.XYZ)
print("Etc, ver documentaciòn")
print(Doblete.S_XYZ)
print(Doblete.T_XYZ)
print(Doblete.OST_XYZ)
print(Doblete.DISTANCE)
print(Doblete.OP)
print(Doblete.TOP)
print(Doblete.TOP_S)
print(Doblete.ALPHA)
print(Doblete.S_LMN)
print(Doblete.LMN)
print(Doblete.R_LMN)
print(Doblete.N0)
print(Doblete.N1)
print(Doblete.WAV)
print(Doblete.G_LMN)
print(Doblete.ORDER)
print(Doblete.GRATING)
print(Doblete.RP)
print(Doblete.RS)
print(Doblete.TP)
print(Doblete.TS)
print(Doblete.TTBE)
print(Doblete.TT)
print(Doblete.BULK_TRANS)
