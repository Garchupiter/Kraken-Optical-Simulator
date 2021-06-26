#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

import time

import numpy as np

import Kraken as kn

start_time = time.time()

config_1 = kn.Kraken_setup()

##############################################################    
P_Obj = kn.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

L1a = kn.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

L1b = kn.surf()
L1b.Rc = (-3.071608670000159E+001)

L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

L1c = kn.surf()
L1c.Rc = (-7.819730726078505E+001)
L1c.Thickness = 9.737604742910693E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

P_Ima = kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 3.0
P_Ima.Name = "Plano imagen"

A = [P_Obj, L1a, L1b, L1c, P_Ima]

#######################################

Doblete = kn.system(A, config_1)

Rayos1 = kn.raykeeper(Doblete)

pSource_0 = [0, 10, 0.0]
dCos = [0.0, 0.0, 1.0]
W = 0.4

Doblete.Trace(pSource_0, dCos, W)
Rayos1.push()
# display2d(Doblete,Rayos1,0)

########################################
print(Doblete.ModificationDetector)
L1a.Rc = (1 + 9.284706570002484E+001)
L1c.Thickness = (1 + 9.737604742910693E+001)
AR, V = Doblete.ModificationDetector()

print(AR)
print(V)
print(np.argwhere(AR != 0.))
