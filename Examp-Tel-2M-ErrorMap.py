# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

import matplotlib.pyplot as plt
import numpy as np

import Kraken as kn

##############################################################   

def ErrorGen():
    L = 1000.
    N = 20.
    hight = 0.001
    SPACE = 2 * L / N
    x = np.arange(-L, L + SPACE, SPACE)
    y = np.arange(-L, L + SPACE, SPACE)
    gx, gy = np.meshgrid(x, y)

    R = np.sqrt((gx * gx) + (gy * gy))

    arg = np.argwhere(R < L)
    Npoints = np.shape(arg)[0]
    X = np.zeros(Npoints)
    Y = np.zeros(Npoints)

    i = 0
    for [a, b] in arg:
        X[i] = gx[a, b]
        Y[i] = gy[a, b]

        i = i + 1

    spa = 10000000
    Z = hight * (np.random.randint(-spa, spa, Npoints)) / (spa * 2.0)

    return [X, Y, Z, SPACE]


##############################################################


P_Obj = kn.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 3500
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1.059E+003 * 2.0

Thickness = 3.452200000000000E+003
M1 = kn.surf()
M1.Rc = -9.638000000004009E+003
M1.Thickness = -Thickness
M1.k = -1.077310000000000E+000
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
M1.InDiameter = 250 * 2.0
M1.Error_map = ErrorGen()

M2 = kn.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness + 1.037525880125084E+003
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0
M2.AxisMove = 0

P_Ima = kn.surf()
P_Ima.Diameter = 1000.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

A = [P_Obj, M1, M2, P_Ima]

######################


configuracion_1 = kn.Kraken_setup()

Telescopio = kn.system(A, configuracion_1)
Rayos = kn.raykeeper(Telescopio)

tam = 9
rad = 2100 / 2

tsis = len(A) - 1

import time

start_time = time.time()

for i in range(-tam, tam + 1):
    for j in range(-tam, tam + 1):

        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            print(i)
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]

            W = 0.4
            Telescopio.Trace(pSource_0, dCos, W)
            Rayos.push()


print("--- %s seconds ---" % (time.time() - start_time))
kn.display3d(Telescopio, Rayos, 2)

print(Telescopio.EFFL)

X, Y, Z, L, M, N = Rayos.pick(-1)

plt.plot(X, Y, 'x')

plt.xlabel('numbers')
plt.ylabel('values')

plt.title('Dot Plot : Red Dots')
plt.axis('square')
plt.show()
