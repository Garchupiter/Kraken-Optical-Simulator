#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Perfect Lens"""

import time
import matplotlib.pyplot as plt
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

start_time = time.time()

# ______________________________________#

P_Obj = Kos.surf(Thickness = 100, Diameter = 30.0)
P_Obj.Name = "o"

L = Kos.surf()
L.Thin_Lens = 50.0
L.Thickness = 50.0
L.Diameter = 30.0

S = Kos.surf(Thickness = 50.0, Diameter = 30.0)
S.Name = "F"

P_Ima = Kos.surf(Diameter = 40.0, Name = "i")

A = [P_Obj, L, S, P_Ima]
config_1 = Kos.Setup()

Lens = Kos.system(A, config_1)
Rays = Kos.raykeeper(Lens)


Surf = 1
W = 0.45
AperVal = L.Diameter
AperType = "EPD"
Pupil = Kos.PupilCalc(Lens, Surf, W, AperType, AperVal)
Pupil.Samp = 1
Pupil.Ptype = "fanx"
Pupil.FieldX = 0.0
Pupil.FieldY = 0.0

# """ FieldType = angle or height """
# Pupil.FieldType = "angle"

# x, y, z, L, M, N = Pupil.Pattern2Field()

# for i in range(0, len(x)):
#     pSource_0 = [x[i], y[i], z[i]]
#     dCos = [L[i], M[i], N[i]]
#     Lens.Trace(pSource_0, dCos, W)
#     Rays.push()

# Kos.display2d(Lens, Rays, 1, arrow=1)






""" FieldType = angle or height """
Pupil.FieldType = "height"

x, y, z, L, M, N = Pupil.Pattern2Field()

for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Lens.Trace(pSource_0, dCos, W)
    Rays.push()

Kos.display2d(Lens, Rays, 1, arrow=1)
