#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp TEl 2M Spyder Spot RMS"""

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
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

P_Obj = Kos.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 1000
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1.059E+003 * 2.0

# ______________________________________#

Spider = Kos.surf()
Spider.Rc = 999999999999.0
Spider.Thickness = 3.452229924716749E+003 + 100.0
Spider.Glass = "AIR"
Spider.Diameter = 1.059E+003 * 2.0

# ______________________________________#

plane1 = pv.Plane(center=[0, 0, 0], direction=[0, 0, 1], i_size=30, j_size=2100, i_resolution=10, j_resolution=10)
plane2 = pv.Plane(center=[0, 0, 0], direction=[0, 0, 1], i_size=2100, j_size=30, i_resolution=10, j_resolution=10)
Baffle1 = pv.Disc(center=[0.0, 0.0, 0.0], inner=0, outer=875 / 2.0, normal=[0, 0, 1], r_res=1, c_res=100)
Baffle2 = Baffle1.merge(plane1)
Baffle3 = Baffle2.merge(plane2)

# ______________________________________#

AAA = pv.MultiBlock()
AAA.append(plane1)
AAA.append(plane2)
AAA.append(Baffle1)

# ______________________________________#

Spider.Mask_Shape = AAA
Spider.Mask_Type = 2
Spider.TiltZ = 0

# ______________________________________#

Thickness = 3.452200000000000E+003
M1 = Kos.surf()
M1.Rc = -9.638000000004009E+003
M1.Thickness = -Thickness
M1.k = -1.077310000000000E+000
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0


# ______________________________________#

M2 = Kos.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness + 1.037535322418897E+003  # 1.037525880125084E+003+1
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0

# ______________________________________#

M2.TiltX = -9.657878504276254E-002
M2.DespY = -2.000000000000000E+000
M2.AxisMove = 0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 100.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, Spider, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()

# ______________________________________#

Telescopio = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Telescopio)

# ______________________________________#

tam = 7
rad = 2200 / 2
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
            Telescopio.Trace(pSource_0, dCos, W)
            Rayos.push()
            W = 0.5
            Telescopio.Trace(pSource_0, dCos, W)
            Rayos.push()
            W = 0.6
            Telescopio.Trace(pSource_0, dCos, W)
            Rayos.push()

# ______________________________________#

Kos.display3d(Telescopio, Rayos, 0)
X, Y, Z, L, M, N = Rayos.pick(-1)

# ______________________________________#

plt.plot(X, Y, 'x')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Spot diagram')
plt.axis('square')
plt.show()
