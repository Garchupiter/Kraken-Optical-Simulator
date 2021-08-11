#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

import time

import matplotlib.pyplot as plt
import numpy as np

import Kraken as kn

start_time = time.time()

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
L1b.Rc = -3.071608670000159E+001

L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

POS_ESP = -40
L1c = kn.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001 + POS_ESP
L1c.Glass = "AIR"
L1c.Diameter = 30

Esp90 = kn.surf()
Esp90.Thickness = POS_ESP
Esp90.Glass = "MIRROR"
Esp90.Diameter = 30.0
Esp90.Name = "Espejo a 90 grados"
Esp90.TiltX = 45.
Esp90.AxisMove = 2.

P_Ima = kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 3.0
P_Ima.Name = "Plano imagen"

A = [P_Obj, L1a, L1b, L1c, Esp90, P_Ima]

######################
config_1 = kn.Kraken_setup()

Doblete = kn.system(A, config_1)

Rayos1 = kn.raykeeper(Doblete)
Rayos2 = kn.raykeeper(Doblete)
Rayos3 = kn.raykeeper(Doblete)
RayosT = kn.raykeeper(Doblete)

tam = 10
rad = 10.0
tsis = len(A) - 1

for j in range(-tam, tam + 1):
    for i in range(-tam, tam + 1):

        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0

            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            W = 0.4

            Doblete.Trace(pSource_0, dCos, W)
            Rayos1.push()
            RayosT.push()

            W=0.5
            Doblete.Trace(pSource_0, dCos,W)
            Rayos2.push()
            RayosT.push()

            W=0.6
            Doblete.Trace(pSource_0, dCos,W)
            Rayos3.push()
            RayosT.push()

kn.display3d(Doblete, RayosT, 2)

# X, Y, Z, L, M, N = RayosT.pick(-1)
# plt.plot(X, Z, 'x')
X, Y, Z, L, M, N = Rayos1.pick(-1)
plt.plot(X, Z, 'x')
X, Y, Z, L, M, N = Rayos2.pick(-1)
plt.plot(X, Z, 'x')
X, Y, Z, L, M, N = Rayos3.pick(-1)
plt.plot(X, Z, 'x')



# axis labeling
plt.xlabel('numbers')
plt.ylabel('values')

# figure name
plt.title('Dot Plot : Red Dots')
plt.axis('square')
plt.show()

