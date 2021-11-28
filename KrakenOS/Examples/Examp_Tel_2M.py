#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M"""

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
M2.AxisMove = 0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 1000.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"
A = [P_Obj, M1, M2, P_Ima]

# ______________________________________#

configuracion_1 = Kos.Setup()
Telescopio = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Telescopio)

# ______________________________________#

W = 0.4
sup = 1
AperVal = 2010
AperType = "EPD"
Pup = Kos.PupilCalc(Telescopio, sup, W, AperType, AperVal)
Pup.Samp = 7
Pup.FieldType = "angle"

# ______________________________________#

Pup.FieldX = 0.0
Pup.Ptype = "hexapolar"
xa, ya, za, La, Ma, Na = Pup.Pattern2Field()

# ______________________________________#

Pup.FieldX = 0.5
Pup.FieldY = 0.5
Pup.Ptype = "square"
xb, yb, zb, Lb, Mb, Nb = Pup.Pattern2Field()

# ______________________________________#

Pup.FieldX = -.5
Pup.FieldY = -.5
Pup.Ptype = "rand"
xc, yc, zc, Lc, Mc, Nc = Pup.Pattern2Field()

# ______________________________________#

for i in range(0, len(xa)):
    pSource_0 = [xa[i], ya[i], za[i]]
    dCos = [La[i], Ma[i], Na[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()
# ______________________________________#

for i in range(0, len(xb)):
    pSource_0 = [xb[i], yb[i], zb[i]]
    dCos = [Lb[i], Mb[i], Nb[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()

# ______________________________________#

for i in range(0, len(xc)):
    pSource_0 = [xc[i], yc[i], zc[i]]
    dCos = [Lc[i], Mc[i], Nc[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()

# ______________________________________#

Kos.display3d(Telescopio, Rayos, 2)
X, Y, Z, L, M, N = Rayos.pick(-1)

# ______________________________________#

plt.plot(X, Y, 'x')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Spot Diagram')
plt.axis('square')
plt.show()
