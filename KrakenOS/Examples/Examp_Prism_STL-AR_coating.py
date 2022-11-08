#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Parabole Mirror Shift"""

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

P_Obj = Kos.surf()
P_Obj.Thickness = 5.0
P_Obj.Diameter = 10
P_Obj.Drawing = 0

P = Kos.surf()
P.Thickness = 10.0
P.Diameter = 5



file = r"prism.stl"
Solid = Kos.surf()
Solid.Diameter = 20
Solid.Solid_3d_stl = file
Solid.Glass= "BK7"
Solid.DespY = -5
Solid.AxisMove = 2
Solid.Thickness = -10
#coating = [[T],[R],[A],[W],[THETA]]

# reflectivity
R = [[0.0, 0.0, 0.0],
     [0.0, 0.0, 0.0]]
# absorption
A = [[0.0, 0.0, 0.0],
     [0.0, 0.0, 0.0]]
# wavelength
W = [0.35, 0.45, 0.55]
# angle
THETA = [0, 45]
# anti reflection coating
Solid.Coating =[R, A, W, THETA]
"""Note: this cannot change the total internal reflection """





P2 = Kos.surf()
P2.Thickness = -10.0
P2.Diameter = 5
P2.Glass = "AIR"

Solid2 = Kos.surf()
Solid2.Diameter = 20
Solid2.Solid_3d_stl = file
Solid2.TiltX = 180
Solid2.Glass= "BK7"
Solid2.DespY = -5
Solid2.AxisMove = 2
Solid2.Thickness = 40

P_Ima = Kos.surf()
P_Ima.Glass = "AIR"
P_Ima.Diameter = 7
P_Ima.Drawing = 1

A = [P_Obj, P, Solid, P2, Solid2, P_Ima]
configuracion_1 = Kos.Setup()
SOLID = Kos.system(A, configuracion_1)
Rays = Kos.raykeeper(SOLID)
SOLID.energy_probability = 0




AA = [P_Obj, P, P_Ima]
SOLID2 = Kos.system(AA, configuracion_1)





W= 0.65
Surf, W, AperVal, AperType = 1, W, P.Diameter, "EPD"
P = Kos.PupilCalc(SOLID2, Surf, W, AperType, AperVal)
P.Samp, P.Ptype, P.FieldX, P.FieldType = 5, "fany", 0, "angle"
x, y, z, L, M, N = P.Pattern2Field()


Kos.NsTraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)

Kos.display2d(SOLID, Rays, 0, arrow = 1)
Kos.display3d(SOLID, Rays, 1)

