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
P_Obj.Thickness = 1
P_Obj.Glass = "AIR"
P_Obj.Diameter = 80
P_Obj.Drawing = 0



# ______________________________________#

file = r"FibInt.stl"


Int = Kos.surf()
Int.Diameter = 80
Int.Solid_3d_stl = file
Int.Thickness = 200
Int.Glass = 1.25






# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0
P_Ima.Thickness = 100.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 1000.0
P_Ima.Drawing = 3




# ______________________________________#

A = [P_Obj,  Int, P_Ima]
configuracion_1 = Kos.Setup()

# ______________________________________#

Fibra = Kos.system(A, configuracion_1)
Rays = Kos.raykeeper(Fibra)
Fibra.Next = 1.1



# Fibra.energy_probability = 0
# Fibra.NsLimit = 30


# ______________________________________#
ang = 45
ARR = np.arange(-ang*0.5, ang*0.5, .5)
for i in ARR:

    tet = i
    pSource_0 = [0 , 0 , 100.0]
    dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
    W = 0.5
    Fibra.NsTrace(pSource_0, dCos, W)
    Rays.push()


# ______________________________________#


Kos.display2d(Fibra, Rays, 0)
# Kos.display3d(Fibra, Rays, 2)
