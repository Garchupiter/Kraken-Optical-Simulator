#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Doublet Lens Zernike"""

import pkg_resources
required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("No instalado")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos


import numpy as np
import matplotlib.pyplot as plt
import time

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 1
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

L1a = Kos.surf()
L1a.Rc = 5.513435044607768E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0

# ______________________________________#

L1b = Kos.surf()
L1b.Rc = -4.408716526030626E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# ______________________________________#

L1c = Kos.surf()
L1c.Rc = -2.246906271406796E+002
L1c.Thickness = 9.737871661422000E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

# ______________________________________#

Z = np.zeros(36)
Z[8] = 0.5
Z[9] = 0.5
Z[10] = 0.5
Z[11] = 0.5
Z[12] = 0.5
Z[13] = 0.5
Z[14] = 0.5
Z[15] = 0.5
L1c.ZNK = Z

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 100.0

# ______________________________________#

A = [P_Obj, L1a, L1b, L1c, P_Ima]
configuracion_1 = Kos.Setup()


# ______________________________________#

Doblete = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Doblete)

# ______________________________________#

tam = 5
rad = 12.0
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

# ______________________________________#

Kos.display3d(Doblete, Rayos, 2)
Kos.display2d(Doblete, Rayos, 1)
