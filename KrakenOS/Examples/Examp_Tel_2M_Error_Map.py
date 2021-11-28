# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M Error Map"""

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
import time


# ______________________________________#

def ErrorGen():
    L = 1000.
    N = 20.
    hight = 0.001
    SPACE = 2 * L / N
    x = np.arange(-L, L + SPACE, SPACE)
    y = np.arange(-L, L + SPACE, SPACE)
    gx, gy = np.meshgrid(x, y)
    R = np.sqrt((gx * gx) + (gy * gy))
    arg = np.argwhere(R < L)
    Npoints = np.shape(arg)[0]
    X = np.zeros(Npoints)
    Y = np.zeros(Npoints)
    i = 0
    for [a, b] in arg:
        X[i] = gx[a, b]
        Y[i] = gy[a, b]
        i = i + 1
    spa = 10000000
    Z = hight * (np.random.randint(-spa, spa, Npoints)) / (spa * 2.0)
    return [X, Y, Z, SPACE]


# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 3500
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
M1.Error_map = ErrorGen()

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

# ______________________________________#

A = [P_Obj, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()

# ______________________________________#

Telescopio = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Telescopio)

# ______________________________________#

tam = 9
rad = 2100 / 2
tsis = len(A) - 1

# ______________________________________#

start_time = time.time()

# ______________________________________#

for i in range(-tam, tam + 1):
    for j in range(-tam, tam + 1):
        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            print(i)
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            W = 0.4
            Telescopio.Trace(pSource_0, dCos, W)
            Rayos.push()

# ______________________________________#

print("--- %s seconds ---" % (time.time() - start_time))
Kos.display3d(Telescopio, Rayos, 2)
print(Telescopio.EFFL)
X, Y, Z, L, M, N = Rayos.pick(-1)

# ______________________________________#

plt.plot(X, Y, 'x')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Spot Diagram')
plt.axis('square')
plt.show()
