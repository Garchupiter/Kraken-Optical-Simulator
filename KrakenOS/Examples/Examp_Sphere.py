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
import os

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Thickness = 1000.0
P_Obj.Diameter = 300
P_Obj.Drawing = 0

# ______________________________________#

M1 = Kos.surf()
M1.Rc = -2 * P_Obj.Thickness
M1.Thickness = M1.Rc / 2
M1.k = -1.0 * 0
M1.Glass = "MIRROR"
M1.Diameter = 2000
M1.CoatingMet = 0



# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Glass = "AIR"
P_Ima.Diameter = 1600.0
P_Ima.Drawing = 0
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, M1, P_Ima]
configuracion_1 = Kos.Setup()

GLASCAT_PATH = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)) + "/KrakenOS/Cat"
MATERIAL_PATH = os.path.join(GLASCAT_PATH, 'Gold.csv')

configuracion_1.LoadMetal(MATERIAL_PATH, "Gold", 1)

# ______________________________________#

Espejo = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Espejo)

# ______________________________________#

W = 0.5876
Surf = 1
AperVal = M1.Diameter
AperType = "EPD"
fieldType = "angle"

Pup = Kos.PupilCalc(Espejo, Surf, W, AperType, AperVal)
Pup.Samp =1
Pup.Ptype = "fany"
Pup.FieldY = 0
x, y, z, L, M, N = Pup.Pattern2Field()
Rayos = Kos.raykeeper(Espejo)




Kos.TraceLoop(x, y, z, L, M, N, W, Rayos, clean = 1)
Kos.display3d(Espejo, Rayos, 0)


# print(Rayos.TP)
# print(Rayos.TS)
print(Rayos.RS[0], Rayos.RS[1], Rayos.RS[2])
print(Rayos.RP[0], Rayos.RP[1], Rayos.RP[2])


# Kos.n_wave_dispersion(configuracion_1, "BK7", 0.4)

# configuracion_1.NAMES
# configuracion_1.CAT
# configuracion_1.NAMES
# configuracion_1.NM
# configuracion_1.ED
# configuracion_1.CD
# configuracion_1.TD
# configuracion_1.OD
# configuracion_1.LD
# configuracion_1.IT



# [Dispersion_Formula, MIL, Nd, Vd, *_] = NM[r]


# x,y,z,l,m,n = Rayos.pick(-1, coordinates="local")







