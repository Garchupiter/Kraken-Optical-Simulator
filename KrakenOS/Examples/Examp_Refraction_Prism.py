#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Diffraction Grating Transmission"""

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

P_Obj = Kos.surf(Thickness = 10, Diameter = 30.0)

Prism_f1 = Kos.surf()
Prism_f1.Rc = 0.0
Prism_f1.Thickness = 13
Prism_f1.Glass = "BK7"
Prism_f1.Diameter = 30.0
Prism_f1.TiltX=20
Prism_f1.AxisMove = 0

Prism_f2 = Kos.surf()
Prism_f2.Rc = 0.0
Prism_f2.Thickness = 10
Prism_f2.Glass = "AIR"
Prism_f2.Diameter = 30.0
Prism_f2.TiltX = -20
Prism_f2.AxisMove = 0

L1_fa = Kos.surf(Rc = 55.1343, Thickness = 6.0, Glass = "BK7", Diameter = 30.0)
L1_fa.DespY = -6.0
L1_fa.TiltX = 25

L1_fb = Kos.surf(Rc = -44.0871, Thickness = 3.0, Glass = "F2", Diameter = 30)
L1_fc = Kos.surf(Rc = -224.690, Thickness = 97.3787, Diameter = 30)

P_Ima = Kos.surf()
P_Ima.Name = "Image plane"
P_Ima.Nm_Pos = [-20,10]
P_Ima.Diameter = 30.0
P_Ima.NumLabel = 0

A = [P_Obj, Prism_f1, Prism_f2, L1_fa, L1_fb, L1_fc, P_Ima]
# A = [P_Obj, Prism_f1, Prism_f2,  P_Ima]
Config = Kos.Setup()

Prism = Kos.system(A, Config)
Rays = Kos.raykeeper(Prism)

# _________________________________________#


n_rays = 10
semidiam = 5.0
nr = int(n_rays/2)
tsis = len(A) - 1
wVs = [0.4, 0.5, 0.6]
for i in range(-nr, nr + 1):
    for j in range(-nr, nr + 1):
        x_0 = (i / nr) * semidiam
        y_0 = (j / nr) * semidiam
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < semidiam:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            for wV in wVs:
                Prism.Trace(pSource_0, dCos, wV)
                Rays.push()

Kos.display2d(Prism, Rays, 0)







