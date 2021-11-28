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

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 50
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

L1a = Kos.surf()
L1a.Thin_Lens = 100.
L1a.Thickness = (100 + 50)
L1a.Rc = 0.0
L1a.Glass = "AIR"
L1a.Diameter = 30.0

# ______________________________________#

L1b = Kos.surf()
L1b.Thin_Lens = 50.
L1b.Thickness = 100.
L1b.Rc = 0.0
L1b.Glass = "AIR"
L1b.Diameter = 30.0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 100.0
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, L1a, L1b, P_Ima]
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
            Doblete.Trace(pSource_0, dCos, W)
            Rayos1.push()
            RayosT.push()

# ______________________________________#

Kos.display3d(Doblete, RayosT, 0)
X, Y, Z, L, M, N = Rayos1.pick(-1)

# ______________________________________#

plt.plot(X, Y, 'x')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Stop Diagram')
plt.axis('square')
plt.show()

# ______________________________________#

print("--- %s seconds ---" % (time.time() - start_time))
