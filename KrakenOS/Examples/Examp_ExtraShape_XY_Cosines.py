#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Extra Shape XY Cosines"""

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

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

L1a = Kos.surf()
L1a.Rc = 55.134 * 0
L1a.Thickness = 9.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0

# ______________________________________#

L1c = Kos.surf()
L1c.Thickness = 40
L1c.Glass = "AIR"
L1c.Diameter = 30


# ______________________________________#

def f(x, y, E):
    r = np.sqrt((x * x) + (y * y * 0))
    H = 2.0 * np.pi * r / E[0]
    zx = np.abs(np.cos(H) * E[1])
    r = np.sqrt((x * x * 0) + (y * y))
    H = 2.0 * np.pi * r / E[0]
    zy = np.abs(np.cos(H) * E[1])
    return zx + zy


# ______________________________________#

coef = [10.0, 1.]
L1c.ExtraData = [f, coef]
L1c.Res = 1

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

Kos.display3d(Lens, Rays, 2)
Kos.display2d(Lens, Rays, 0)
