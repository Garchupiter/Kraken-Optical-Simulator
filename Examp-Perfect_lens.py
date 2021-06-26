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

config_1 = kn.Kraken_setup()

##############################################################


P_Obj = kn.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 50
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

L1a = kn.surf()
L1a.Thin_Lens = 100.
L1a.Thickness = (100 + 50)
L1a.Rc = 0.0

L1a.Glass = "AIR"
L1a.Diameter = 30.0

L1b = kn.surf()
L1b.Thin_Lens = 50.
L1b.Thickness = 100.
L1b.Rc = 0.0

L1b.Glass = "AIR"
L1b.Diameter = 30.0

P_Ima = kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 100.0
P_Ima.Name = "Plano imagen"

A = [P_Obj, L1a, L1b, P_Ima]

######################


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

            # W=0.5
            # Doblete.Trace(pSource_0, dCos,W)
            # Rayos2.push()
            # RayosT.push()

            # W=0.6
            # Doblete.Trace(pSource_0, dCos,W)
            # Rayos3.push()
            # RayosT.push()

kn.display3d(Doblete, RayosT, 0)

# #display2d(Doblete,Rayos1,0)

# print(Doblete.EFFL)


X, Y, Z, L, M, N = Rayos1.pick(-1)
plt.plot(X, Y, 'x')

# plotpoints1=Rayos1(Rayos2)
# X,Y,Z,L,M,N=plotpoints1.pick(-1)
# plt.plot(X,Y, 'x')

# plotpoints2=Rayos1(Rayos3)
# X,Y,Z,L,M,N=plotpoints2.pick(-1)
# plt.plot(X,Y, 'x')


# axis labeling
plt.xlabel('numbers')
plt.ylabel('values')

# figure name
plt.title('Dot Plot : Red Dots')
plt.axis('square')
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
