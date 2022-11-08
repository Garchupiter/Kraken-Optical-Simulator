#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Flat Mirror 45 Deg"""

import time
import matplotlib.pyplot as plt
import numpy as np

import sys
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
# ______________________________________#

start_time = time.time()

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

L1a = Kos.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# ______________________________________#

L1b = Kos.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# ______________________________________#

POS_ESP = -40
L1c = Kos.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001 + POS_ESP
L1c.Glass = "AIR"
L1c.Diameter = 30

# ______________________________________#

Esp90 = Kos.surf()
Esp90.Thickness = POS_ESP
Esp90.Glass = "MIRROR"
Esp90.Diameter = 30.0
Esp90.Name = "Espejo a 90 grados"
Esp90.TiltX = 45.
Esp90.AxisMove = 2.



# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 3.0
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, L1a, L1b, L1c, Esp90, P_Ima]

config_1 = Kos.Setup()

# ______________________________________#

Doblete = Kos.system(A, config_1)


Rayos1 = Kos.raykeeper(Doblete)
Rayos2 = Kos.raykeeper(Doblete)
Rayos3 = Kos.raykeeper(Doblete)
RayosT = Kos.raykeeper(Doblete)

# ______________________________________#

tam = 10
rad = 10.0
tsis = len(A) - 1
for j in range(-tam, tam + 1):
    for i in range(-tam, tam + 1):
        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            W = 0.4
            Doblete.FastTrace(pSource_0, dCos, W)
            Rayos1.push()
            RayosT.push()
            W = 0.5
            Doblete.FastTrace(pSource_0, dCos, W)
            Rayos2.push()
            RayosT.push()
            W = 0.6
            Doblete.FastTrace(pSource_0, dCos, W)
            Rayos3.push()
            RayosT.push()

# ______________________________________#

# Kos.display3d(Doblete, RayosT, 2)

system = [Doblete, Doblete2]
Kos.display2dPlus(system, RayosT, 0)


# # ______________________________________#

# X, Y, Z, L, M, N = Rayos1.pick(-1)
# plt.plot(X, Z, 'x')
# X, Y, Z, L, M, N = Rayos2.pick(-1)
# plt.plot(X, Z, 'x')
# X, Y, Z, L, M, N = Rayos3.pick(-1)
# plt.plot(X, Z, 'x')
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Spot Diagram')
# plt.axis('square')
# plt.show()
