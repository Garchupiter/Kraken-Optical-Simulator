#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Extra Shape Micro Lens Array"""

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
import numpy as np
import matplotlib.pyplot as plt

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

L1a = Kos.surf()
L1a.Rc = 55.134 * 0
L1a.Thickness = 2.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0

# ______________________________________#

L1c = Kos.surf()
L1c.Thickness = 40
L1c.Glass = "AIR"
L1c.Diameter = 30


# ______________________________________#

def f(x, y, E):
    DeltaX = E[0] * np.rint(x / E[0])
    DeltaY = E[0] * np.rint(y / E[0])
    x = x - DeltaX
    y = y - DeltaY
    s = np.sqrt((x * x) + (y * y))
    c = 1.0 / E[1]
    InRoot = 1 - (E[2] + 1.0) * c * c * s * s
    z = (c * s * s / (1.0 + np.sqrt(InRoot)))
    return z


# ______________________________________#

coef = [3.0, -3, 0]
L1c.ExtraData = [f, coef]
L1c.Res = 2

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 300.0
P_Ima.Name = "Image plane"

# ______________________________________#

A = [P_Obj, L1a, L1c, P_Ima]
Config_1 = Kos.Setup()

# ______________________________________#

Lens = Kos.system(A, Config_1)
Rays = Kos.raykeeper(Lens)

# ______________________________________#

Wav = 0.45
for i in range(-100, 100 + 1):
    pSource = [0.0, i / 10., 0.0]
    dCos = [0.0, 0.0, 1.0]
    Lens.Trace(pSource, dCos, Wav)
    Rays.push()

# ______________________________________#

Kos.display3d(Lens, Rays, 1)
# Kos.display2d(Lens, Rays, 0)
