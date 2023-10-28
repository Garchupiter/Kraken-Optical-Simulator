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
P_Obj.Rc = 0.0
P_Obj.Thickness = 51
P_Obj.Glass = "AIR"
P_Obj.Diameter = 50.0

# ______________________________________#

f1 = Kos.surf()
f1.Glass = 1.343
f1.TiltX=90
f1.AxisMove = 1
f1.Diameter = 100.0
f1.Thickness=1

f2 = Kos.surf()
f2.Glass = 1.522
f2.Diameter = 100.0
f2.Thickness=1

f3 = Kos.surf()
f3.Glass = 1.343
f3.Diameter = 100.0
f3.Thickness=1

f4 = Kos.surf()
f4.Glass = "AIR"
f4.Diameter = 100.0
f4.Thickness=0


# ______________________________________#
P_Ima = Kos.surf()
P_Ima.Rc = 0
P_Ima.Thickness = 0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 50.0
P_Ima.Drawing = 1

P_Ima.TiltX=-90
P_Ima.DespY =100
P_Ima.AxisMove = 1
P_Ima.Diameter = 100.0

# ______________________________________#

A = [P_Obj, f1, f2, f3, f4, P_Ima]
configuracion_1 = Kos.Setup()
Tube = Kos.system(A, configuracion_1)
RayosT = Kos.raykeeper(Tube)

# ______________________________________#
ang = 10
ARR = np.arange(-ang*0.5, ang*0.5, 0.1)
for i in ARR:

    tet = i
    pSource_0 = [0 , -1.5 , 0.0]
    dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
    W = 0.5
    Tube.NsTrace(pSource_0, dCos, W)
    RayosT.push()


# ______________________________________#


Kos.display2d(Tube, RayosT, 0)
Kos.display3d(Tube, RayosT, 0)
