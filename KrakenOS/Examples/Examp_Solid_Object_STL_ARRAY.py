#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Solid Objects STL ARRAY"""

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
P_Obj.Thickness = 5000.0
P_Obj.Glass = "AIR"
P_Obj.Diameter = 6.796727741707513E+002 * 2.0
P_Obj.Drawing = 0

# ______________________________________#

FOV = 0.5
Ref_esp = 0.8
Tf = 200  # 40
A1 = Tf * Tf
Conc = 800
Conc = Conc / 0.8
fpl = int(np.round(np.sqrt(Conc) / 2.0))
print("Numero de facetas por lado (Calculo 1)", (fpl * 2.0) + 1.0)
print("Tamanio por lado(mm) ", Tf * ((fpl * 2.0) + 1.0))

# ______________________________________#

n = fpl
FN = 1
focal = Tf * (fpl * 2.0 + 1.0) * FN
print("Distancia focal(mm) ", focal)

# ______________________________________#

sobredim = 2.0 * focal * np.tan(np.deg2rad(FOV / 2.0))
print("Sobredim (mm): ", sobredim)

# ______________________________________#

Tf = Tf - sobredim
A2 = Tf * Tf
RA = A1 / A2
Conc = Conc * RA
print("Nuevo tamaño de las facetas: ", Tf)

# ______________________________________#

fpl = int(np.round(np.sqrt(Conc) / 2.0))
print("Nuevo numero de facetas por lado (Calculo 2)", (fpl * 2.0) + 1.0)
print("Tamanio por lado 2a (mm) ", Tf * ((fpl * 2.0) + 1.0))

# ______________________________________#

n = fpl
Cx = Tf
Cy = Tf
Cz = 0
Lx = Tf
Ly = Tf
Lz = 1.0

# ______________________________________#

element0 = pv.Cube(center=(0.0, 0.0, 0.0), x_length=0.1, y_length=0.1, z_length=0.1, bounds=None)
for A in range(-n, n + 1):
    for B in range(-n, n + 1):
        Ty = 0.5 * np.rad2deg(np.arctan2(Cx * A, focal))
        Tx = -0.5 * np.rad2deg(np.arctan2(Cy * B, focal))
        element1 = pv.Cube(center=(0.0, 0.0, 0.0), x_length=Lx, y_length=Ly, z_length=Lz, bounds=None)
        element1.rotate_x(Tx)
        element1.rotate_y(Ty)
        v = [-Cx / 2.0, -Cy / 2.0, 0]
        element1.translate(v)
        v = [Cx * A, Cy * B, Cz]
        element1.translate(v)
        element0 = element0 + element1
element0.save("salida.stl")
direc = r"salida.stl"

# ______________________________________#

objeto = Kos.surf()
objeto.Diameter = 118.0 * 2.0
objeto.Solid_3d_stl = direc
objeto.Thickness = -6000
objeto.Glass = "MIRROR"
objeto.TiltX = 0
objeto.TiltY = 0
objeto.DespX = 0
objeto.DespY = 0
objeto.AxisMove = 0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0
P_Ima.Thickness = -1.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 2000.0
P_Ima.Drawing = 1
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, objeto, P_Ima]
configur = Kos.Setup()

# ______________________________________#

Telescope = Kos.system(A, configur)
Rays = Kos.raykeeper(Telescope)

# ______________________________________#

W = 0.633
tam = 25
rad = 5500.0
tsis = len(A) + 2
for j in range(-tam, tam + 1):
    for i in range(-tam, tam + 1):
        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            Telescope.NsTrace(pSource_0, dCos, W)
            if np.shape(Telescope.NAME)[0] != 0:
                if Telescope.NAME[-1] == "Plano imagen":
                    plt.plot(Telescope.Hit_x[-1], Telescope.Hit_y[-1], '.', c="g")
                    Rays.push()

# ______________________________________#

plt.axis('square')
plt.show()
Kos.display3d(Telescope, Rays, 0)
