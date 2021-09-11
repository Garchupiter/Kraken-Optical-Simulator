#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M Spyder Spot Diagram"""

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import KrakenOS as Kos
import scipy

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

plane1 = pv.Plane(center=[0, 0, 0], direction=[0, 0, 1], i_size=10, j_size=2100, i_resolution=10, j_resolution=10)
plane2 = pv.Plane(center=[0, 0, 0], direction=[0, 0, 1], i_size=2100, j_size=10, i_resolution=10, j_resolution=10)
Baffle1 = pv.Disc(center=[0.0, 0.0, 0.0], inner=0, outer=875 / 2.0, normal=[0, 0, 1], r_res=1, c_res=100)
# Baffle2 = Baffle1.boolean_add(plane1)
# Baffle3 = Baffle2.boolean_add(plane2)
AAA = pv.MultiBlock()
AAA.append(plane1)
AAA.append(plane2)
AAA.append(Baffle1)


Spider.Mask_Shape = AAA
Spider.Mask_Type = 2
Spider.TiltZ = 0

# ______________________________________#

Thickness = 3.452229924716749E+003
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
focusShift = 1.0  # Set cero to focus, 1 is only to see the spider is the spot diagram
M2.Thickness = Thickness + 1.037179115116706E+003 + focusShift
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 100.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, Spider, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()

# ______________________________________#

Telescope = Kos.system(A, configuracion_1)
Rays = Kos.raykeeper(Telescope)

# ______________________________________#

# Gaussian (Seeing)
def f(x):
    x = np.rad2deg(x)
    seing = 1.2 / 3600.0
    sigma = seing / 2.3548
    mean = 0
    standard_deviation = sigma
    y = scipy.stats.norm(mean, standard_deviation)
    res = y.pdf(x)
    return res

Sun = Kos.SourceRnd()
Sun.field = 4 * 1.2 / (2.0 * 3600.0)

Sun.fun = f
Sun.dim = 2100/2.0 # radio del espejo primario
Sun.num = 1000
Sun.type = 0 # 1 for Square, 0 for circle

W=0.7

L, M, N, X, Y, Z = Sun.rays()

Xr = np.zeros_like(L)
Yr = np.zeros_like(L)
Nr = np.zeros_like(L)
Energy = np.zeros_like(L)


con = 0
con2 = 0

for i in range(0, Sun.num):
    if con2 == 20:
        print(100. * i / Sun.num)
        con2 = 0

    pSource_0 = [X[i], Y[i], Z[i]]
    dCos = [L[i], M[i], N[i]]
    Telescope.Trace(pSource_0, dCos, W)

    Xr[con] = Telescope.Hit_x[-1]
    Yr[con] = Telescope.Hit_y[-1]

    Nr[con] = Telescope.SURFACE[-1]

    # print(Telescope.TT)
    Energy[con] = Telescope.TT

    con = con + 1
    con2 = con2 + 1
    # Rays.push()

args = np.argwhere(Nr == len(A)-1)

plt.plot(Xr[args], Yr[args], '.', c="g", markersize=1)

# axis labeling
plt.xlabel('x')
plt.ylabel('y')

# figure name
plt.title('Spot diagram')
plt.axis('square')
plt.show()

#             Rays.push()
# Kos.display3d(Telescope, Rays, 0)

print("----------------------------")
print("Numero de rayos trasados:")
print(len(Nr))

print("Numero de rayos que llegan al plano imagen:")
print(len(args))

print("Energia total en rayos que llegan al plano imagen:")
En=Energy[args]
EnT=np.sum(En)
print(EnT)

