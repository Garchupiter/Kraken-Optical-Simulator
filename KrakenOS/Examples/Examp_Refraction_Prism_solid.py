#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Diffraction Grating Transmission"""

import numpy as np
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
import pyvista as pv




# mesh points
r = 10
hx = r/2
rad = np.pi / 180

p0x = r/2
p0y = r * np.sin(((90 + 120)*rad))
p0z = r * np.cos(((90 + 120)*rad))

p1x = r/2
p1y = r * np.sin((90 * rad))
p1z = r * np.cos((90 * rad))

p2x = r/2
p2y = r * np.sin(((90 - 120)*rad))
p2z = r * np.cos(((90 - 120)*rad))

p3x = -p0x
p3y = p0y
p3z = p0z

p4x = -p1x
p4y = p1y
p4z = p1z

p5x = -p2x
p5y = p2y
p5z = p2z


H = -r/2
vertices = np.array([

   [p0x, p0y + H, p0z], [p1x, p1y + H, p1z], [p2x, p2y + H, p2z],
   [p3x, p3y + H, p3z], [p4x, p4y + H, p4z], [p5x, p5y + H, p5z]
                  ])

# mesh faces
faces = np.hstack(
    [
        [3, 0, 1, 2],  # triangle
        [3, 3, 4, 5],  # triangle
        [4, 0, 1, 4, 3],  # rectangle
        [4, 1, 2, 5, 4],  # rectangle
        [4, 0, 2, 5, 3],  # rectangle

    ]
)

Solid = pv.PolyData(vertices, faces)


P_Obj = Kos.surf(Thickness = 30, Diameter = 5.0)

Prism_SolObj = Kos.surf()
Prism_SolObj.Thickness = 30
Prism_SolObj.Glass = "BK7"
Prism_SolObj.Solid_3d_stl = Solid
Prism_SolObj.TiltX = 16
Prism_SolObj.AxisMove = 1.99

R = [[0.0, 0.0, 0.0],
     [0.0, 0.0, 0.0]]

A = [[0.0, 0.0, 0.0],
     [0.0, 0.0, 0.0]]

W = [0.4, 0.5, 0.6]

THETA = [0, 45]
Prism_SolObj.Coating =[R, A, W, THETA]

P_Ima = Kos.surf(Diameter = 20.0, NumLabel = 0)

A = [P_Obj, Prism_SolObj, P_Ima]
Config = Kos.Setup()

Prism = Kos.system(A, Config)
Rays = Kos.raykeeper(Prism)

# _________________________________________#

Prism.energy_probability = 0
n_rays = 5
semidiam = 0.1
nr = int(n_rays/2)
tsis = len(A) - 1
wVs = [0.4, 0.5, 0.6]
for i in range(-nr, nr + 1):
    for j in range(-nr, nr + 1):
        x_0 = (i / nr) * semidiam
        y_0 = (j / nr) * semidiam
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < semidiam:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            for wV in wVs:
                Prism.NsTrace(pSource_0, dCos, wV)
                Rays.push()

Kos.display2d(Prism, Rays, 0)







