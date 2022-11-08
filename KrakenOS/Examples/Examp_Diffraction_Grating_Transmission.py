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
# _________________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# _________________________________________#

Dif_Obj_c1 = Kos.surf()
Dif_Obj_c1.Rc = 0.0
Dif_Obj_c1.Thickness = 1
Dif_Obj_c1.Glass = "BK7"
Dif_Obj_c1.Diameter = 30.0
# Dif_Obj_c1.Grating_D = 1.0
# Dif_Obj_c1.Diff_Ord = 1.
# Dif_Obj_c1.Grating_Angle = 45.

# _________________________________________#

Dif_Obj_c2 = Kos.surf()
Dif_Obj_c2.Rc = 0.0
Dif_Obj_c2.Thickness = 10
Dif_Obj_c2.Glass = "AIR"
Dif_Obj_c2.Diameter = 30.0

# _________________________________________#

L1a = Kos.surf()
L1a.Rc = 5.513435044607768E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0

# _________________________________________#

L1b = Kos.surf()
L1b.Rc = -4.408716526030626E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# _________________________________________#

L1c = Kos.surf()
L1c.Rc = -2.246906271406796E+002
L1c.Thickness = 9.737871661422000E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

# _________________________________________#

P_Ima = Kos.surf()
P_Ima.Name = "Plano imagen"
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 300.0

# _________________________________________#

A = [P_Obj, Dif_Obj_c1, Dif_Obj_c2, L1a, L1b, L1c, P_Ima]
configuracion_1 = Kos.Setup()

# _________________________________________#

Doblete = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Doblete)

# _________________________________________#

tam = 5
rad = 10.0
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
            Doblete.Trace(pSource_0, dCos, W)
            Rayos.push()
            W = 0.5
            Doblete.Trace(pSource_0, dCos, W)
            Rayos.push()
            W = 0.6
            Doblete.Trace(pSource_0, dCos, W)
            Rayos.push()

# _________________________________________#

Kos.display3d(Doblete, Rayos, 0)
