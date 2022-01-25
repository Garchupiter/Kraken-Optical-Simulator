
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

from scipy import optimize
from scipy.optimize import fmin_cg

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

import numpy as np
import matplotlib.pyplot as plt
import time


P_Obj=Kos.surf()
P_Obj.Rc=0
P_Obj.Thickness=3000+3.452200000000000E+003
P_Obj.Glass="AIR"
P_Obj.Diameter=6502.4
P_Obj.Drawing=0

Thickness=6178
M1=Kos.surf()
M1.Rc=-16256.0
M1.Thickness=-Thickness
M1.k=-1.
M1.Glass="MIRROR"
M1.Diameter=6502.4

M2=Kos.surf()
M2.Rc=-5150.974
M2.Thickness=Thickness+1851.28
M2.k=-2.6946
M2.Glass="MIRROR"

M2.Diameter=1714.5
M2.DespX=0
M2.AxisMove=0

P_Ima=Kos.surf()
P_Ima.Diameter=1000.0
P_Ima.Glass="AIR"
P_Ima.Name="Plano imagen"

A = [P_Obj, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()
Telescope = Kos.system(A, configuracion_1)

# ______________________________________#

W=0.4861
pSource_0 = [6500/2., 0, 0]
dCos = [0., 0., 1.]
Telescope.Trace(pSource_0, dCos, W)
print(1.0)
for i in range(1,4):
    P=Telescope.OST_XYZ[-i]
    print(P[0])

print(" ")
for i in range(1,4):
    P=Telescope.OST_XYZ[-i]
    print(P[1])


print(type(P[1]))
