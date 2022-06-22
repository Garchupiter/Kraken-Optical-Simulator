#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens Tilt nonSec"""

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

L1a = Kos.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.AxisMove = 1
L1a.TiltX = 13.0
L1a.DespZ = 5.0

# reflectivity
R = [[0.5, 0.5, 0.5],
     [0.0, 0.0, 0.0]]
# absorption
A = [[0.0, 0.0, 0.0],
     [0.0, 0.0, 0.0]]
# wavelength
W = [0.35, 0.45, 0.55]
# angle
THETA = [0, 45]
# anti reflection coating
L1a.Coating =[R, A, W, THETA]
"""Note: this cannot change the total internal reflection """

# _________________________________________#

L1b = Kos.surf()
L1b.Rc = (-3.071608670000159E+001)
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30


THETA = [0, 45]

L1b.Coating =[R, A, W, THETA]

# _________________________________________#

L1c = Kos.surf()
L1c.Rc = (-7.819730726078505E+001)
L1c.Thickness = 9.737604742910693E+001
L1c.Glass = "AIR"
L1c.Diameter = 30


L1c.Coating =[R, A, W, THETA]


# _________________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 1000.0
P_Ima.Name = "Plano imagen"

# _________________________________________#

A = [P_Obj, L1a, L1b, L1c, P_Obj, L1a, L1b, L1c, P_Ima]
configuracion_1 = Kos.Setup()

# _________________________________________#

Doblete = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Doblete)
Doblete.energy_probability=1

# _________________________________________#

tam = 30
rad = 18.0
tsis = len(A) - 1
for j in range(-tam, tam + 1):
    i = 0
    x_0 = (i / tam) * rad
    y_0 = (j / tam) * rad
    r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
    if r < rad:
        tet = 0.0
        pSource_0 = [x_0, y_0, 0.0]
        dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
        W = 0.4
        Doblete.NsTrace(pSource_0, dCos, W)
        Rayos.push()
        W = 0.5
        Doblete.NsTrace(pSource_0, dCos, W)
        Rayos.push()
        W = 0.6
        Doblete.NsTrace(pSource_0, dCos, W)
        Rayos.push()

# _________________________________________#

Kos.display3d(Doblete, Rayos, 2)
