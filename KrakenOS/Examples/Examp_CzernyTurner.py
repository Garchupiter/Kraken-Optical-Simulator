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


F = 300
D = 50
h = 50
k = -1

r = np.sqrt((F**2.0) + (h**2.0))
Theta=np.rad2deg(np.arctan(h/F))
nh = r*(h/F)

W =0.65
d=1
m=1
Tm = np.arcsin((m*W/(2.0*d)))
Tm = np.rad2deg(Tm)


P_Obj = Kos.surf()
P_Obj.Thickness = r
P_Obj.Diameter = D
P_Obj.Drawing = 0

OAP1 = Kos.surf()
OAP1.Rc = -2.0 * F
OAP1.k  = k
OAP1.Thickness = -r
OAP1.Glass = "MIRROR"
OAP1.Diameter = D
OAP1.ShiftY = -h
OAP1.TiltX = -Theta
OAP1.Order = 0.0
OAP1.AxisMove = 0.0

M = Kos.surf()
M.Thickness = r
M.Glass = "MIRROR"
M.Diameter = D
M.DespY = -nh
M.TiltX = Tm
M.Order = 0.0
M.AxisMove = 0.0
M.Grating_D = d
M.Diff_Ord = m
M.Grating_Angle = 0.0

OAP2 = Kos.surf()
OAP2.Rc = -2.0 * F
OAP2.k  = k
OAP2.Thickness = - r
OAP2.Glass = "MIRROR"
OAP2.Diameter = D
OAP2.ShiftY = h
OAP2.DespY = -nh * 2.0
OAP2.TiltX = Theta
OAP2.Order = 0.0
OAP2.AxisMove = 0.0

P_Ima = Kos.surf()
P_Ima.DespY = -nh * 2.0
P_Ima.TiltX = -Theta/2.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 20

A = [P_Obj, OAP1, M, OAP2, P_Ima]
configuracion_1 = Kos.Setup()

ZernTurn = Kos.system(A, configuracion_1)
Rays = Kos.raykeeper(ZernTurn)

Surf, W, AperVal, AperType = 1, W, OAP1.Diameter, "EPD"
P = Kos.PupilCalc(ZernTurn, Surf, W, AperType, AperVal)
P.Samp, P.Ptype, P.FieldX, P.FieldType = 4, "fany", 0, "height"
x, y, z, L, M, N = P.Pattern2Field()

delta=0.001
Rang= 0.01
DW = np.arange(-Rang,Rang+delta,delta)
for dw in DW:
    Kos.TraceLoop(x, y, z, L, M, N, W + dw, Rays, clean = 0)




W =0.35
Rays2 = Kos.raykeeper(ZernTurn)

Surf, W, AperVal, AperType = 1, W, OAP1.Diameter, "EPD"
P = Kos.PupilCalc(ZernTurn, Surf, W, AperType, AperVal)
P.Samp, P.Ptype, P.FieldX, P.FieldType = 4, "fany", 0, "height"
x, y, z, L, M, N = P.Pattern2Field()

delta=0.001
Rang= 0.01
DW = np.arange(-Rang,Rang+delta,delta)
for dw in DW:
    Kos.TraceLoop(x, y, z, L, M, N, W + dw, Rays2, clean = 0)




# Kos.display3d(ZernTurn, Rays, 0)


P3D = Kos.display3d_OB()

P3D.SYSTEM = ZernTurn
P3D.RAYS = Rays
P3D.plot()



P3D.RAYS = Rays2
P3D.plot()


# d=1
# Tii= Theta
# Ti = np.deg2rad(Tii)
# m=1
# Tm = np.arcsin((m*W/(2.0*d)))

# Tm = np.rad2deg(Tm)
# print(Tm)
