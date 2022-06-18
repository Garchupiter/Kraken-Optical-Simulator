#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M Spyder Spot Diagram"""

# import os
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

P_Obj = Kos.surf()
P_Obj.Thickness = 2.000000000000000E+003
P_Obj.Glass = "AIR"
P_Obj.Diameter = 6.796727741707513E+002 * 2.0
P_Obj.Drawing = 0

# ______________________________________#

M1 = Kos.surf()
M1.Rc = -6.0E+003
M1.Thickness = -1.774190000000000E+003 + 1.853722901194000E+000
M1.k = -1.6+000
M1.Glass = "MIRROR"
M1.Diameter = 6.63448E+002 * 2.0
M1.InDiameter = 228.6 * 2.0
M1.DespY = 0.0
M1.TiltX = 0.0000
M1.AxisMove = 1

# ______________________________________#

M2 = Kos.surf()
M2.Rc = -6.00E+003
M2.Thickness = -M1.Thickness
M2.k = -3.4782E+001
M2.Glass = "MIRROR"
M2.Diameter = 3.0E+002 * 2.0

# ______________________________________#

Vertex = Kos.surf()
Vertex.Thickness = 130.0
Vertex.Glass = "AIR"
Vertex.Diameter = 600.0
Vertex.Drawing = 0

# ______________________________________#

file = r"Prisma.stl"

# ______________________________________#

objeto = Kos.surf()
objeto.Diameter = 118.0 * 2.0
objeto.Solid_3d_stl = file
objeto.Thickness = 600
objeto.Glass = "BK7"
objeto.TiltX = 55
objeto.TiltY = 0
objeto.TiltZ = 45
objeto.DespX = 0
objeto.DespY = 0
objeto.AxisMove = 0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0
P_Ima.Thickness = 100.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 500.0
P_Ima.Drawing = 1

# ______________________________________#

A = [P_Obj, M1, M2, Vertex, objeto, P_Ima]
configuracion_1 = Kos.Setup()

# ______________________________________#

Telescope = Kos.system(A, configuracion_1)
Rays = Kos.raykeeper(Telescope)




Telescope.energy_probability = 1
Telescope.NsLimit = 10

# ______________________________________#

W = 0.633
tam = 5
rad = 6.56727741707513E+002
tsis = len(A) + 2
for gg in range(0, 10):
    for j in range(-tam, tam + 1):
        # j=0
        for i in range(-tam, tam + 1):
            x_0 = (i / tam) * rad
            y_0 = (j / tam) * rad
            r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
            if r < rad:
                tet = 0.0
                pSource_0 = [x_0, y_0, 0.0]
                # print("-...............")
                dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
                W = 0.633
                Telescope.NsTrace(pSource_0, dCos, W)
                Rays.push()

# # ______________________________________#

Kos.display3d(Telescope, Rays, 0)
print(Telescope.EFFL)
