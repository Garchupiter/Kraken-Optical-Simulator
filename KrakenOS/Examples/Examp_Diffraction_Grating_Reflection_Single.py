#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Diffraction Grating Reflection"""

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
# _________________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 250
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0


# _________________________________________#

Dif_Obj = Kos.surf()
Dif_Obj.Rc = 0.0
Dif_Obj.Thickness = -250
Dif_Obj.Glass = "MIRROR"
Dif_Obj.Diameter = 30.0
Dif_Obj.Grating_D = 1000/600
Dif_Obj.Diff_Ord = 1
Dif_Obj.Grating_Angle = 90.0


# _________________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Name = "Plano imagen"
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 1000.0
P_Ima.Drawing = 0

# _________________________________________#

A = [P_Obj, Dif_Obj, P_Ima]
configuracion_1 = Kos.Setup()

# _________________________________________#

Doblete = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Doblete)

# _________________________________________#

pSource_0 = [0, 0, 0.0]
tet = 0
dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
W = 0.650
Doblete.Trace(pSource_0, dCos, W)
print(Doblete.XYZ[-1])
Rayos.push()


# ______________________________________#

Kos.display3d(Doblete, Rayos, 1)
