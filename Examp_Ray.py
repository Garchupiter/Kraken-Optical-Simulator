#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Ray"""

import numpy as np
import KrakenOS as Kn

# ______________________________________#

P_Obj = Kn.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 0.1
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

P_Obj2 = Kn.surf()
P_Obj2.Rc = 0.0
P_Obj2.Thickness = 10
P_Obj2.Glass = "AIR"
P_Obj2.Diameter = 30.0

# ______________________________________#

L1a = Kn.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# ______________________________________#

L1b = Kn.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# ______________________________________#

L1c = Kn.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

# ______________________________________#

P_Ima = Kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 18.0
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, P_Obj2, L1a, L1b, L1c, P_Ima]
configuracion_1 = Kn.Kraken_setup()

# ______________________________________#

Doblete = Kn.system(A, configuracion_1)
Rayos = Kn.raykeeper(Doblete)

# ______________________________________#

pSource_0 = [0, 14, 0]
tet = 0.1
dCos = [0.0, np.sin(np.deg2rad(tet)), -np.cos(np.deg2rad(tet))]
W = 0.4
Doblete.Trace(pSource_0, dCos, W)
Rayos.push()

# ______________________________________#

Kn.display3d(Doblete, Rayos, 2)
