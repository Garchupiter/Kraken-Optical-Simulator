# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M Wavefront Fitting"""

import os
import sys
import matplotlib.pyplot as plt
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

# ______________________________________#

currentDirectory = os.getcwd()
sys.path.insert(1, currentDirectory + '/library')

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 1000
P_Obj.Glass = "AIR"
P_Obj.Diameter = 150.0 * 2.0
P_Obj.Drawing=1

# ______________________________________#


N0 = Kos.surf()
N0.Glass = "AIR"
N0.Diameter = 300
N0.Drawing=0

#############################
N1 = Kos.surf()
N1.Glass = "NULL"
N1.TiltZ = 0
N1.AxisMove = 1
N1.Drawing=0

#############################

N2 = Kos.surf()
N2.Glass = "NULL"
N2.TiltY = 20.7
N2.AxisMove = 1
N2.Drawing=0

###########################

G1 = Kos.surf()
G1.Thickness = -1000
G1.Glass = "MIRROR"
G1.Diameter = 1300
G1.Grating_D =1/0.3
G1.Diff_Ord = -1
G1.Grating_Angle = 90


P_Ima = Kos.surf()
P_Ima.Diameter = 498*4
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj,N0, N1, N2, G1, P_Ima]
configuracion_1 = Kos.Setup()
Telescopio = Kos.system(A, configuracion_1)

# ______________________________________#
Telescopio = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Telescopio)

# ______________________________________#

tam = 9
rad = 100 / 2
tsis = len(A) - 1

# ______________________________________#



# ______________________________________#




# tet = 0.0
# pSource_0 = [0, 0, 0.0]
# dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
# W = .55
# Telescopio.Trace(pSource_0, dCos, W)

# # print(Telescopio.Hit_x[-1], Telescopio.Hit_y[-1])

# print(Telescopio.OST_XYZ[-1])

# Rayos.push()






# for i in range(-tam, tam + 1):
for j in range(-tam, tam + 1):
    x_0 = (0 / tam) * rad
    y_0 = (j / tam) * rad
    r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
    if r < rad:
        print(0)
        tet = 4.28*2
        pSource_0 = [x_0, y_0, 0.0]
        dCos = [0.0, -np.sin(np.deg2rad(tet)), -np.cos(np.deg2rad(tet))]
        print(dCos)
        print([-0.00314876, -0.07473735, -0.99719828])
        # dCos = [-0.00314876, -0.07473735, -0.99719828]
        W = .55
        Telescopio.Trace(pSource_0, dCos, W)
        Rayos.push()



for j in range(-tam, tam + 1):
    x_0 = (0 / tam) * rad
    y_0 = (j / tam) * rad
    r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
    if r < rad:
        print(0)
        tet = -4.28*2
        pSource_0 = [x_0, y_0, 0.0]
        dCos = [0.0, -np.sin(np.deg2rad(tet)), -np.cos(np.deg2rad(tet))]
        print(dCos)
        print([-0.00314876, -0.07473735, -0.99719828])
        # dCos = [-0.00314876, -0.07473735, -0.99719828]
        W = .55
        Telescopio.Trace(pSource_0, dCos, W)
        Rayos.push()







for j in range(-tam, tam + 1):
    x_0 = (0 / tam) * rad
    y_0 = (j / tam) * rad
    r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
    if r < rad:
        print(0)
        tet = 0.0
        pSource_0 = [x_0, y_0, 0.0]
        dCos = [0.0, -np.sin(np.deg2rad(tet)), -np.cos(np.deg2rad(tet))]
        print(dCos)
        print([-0.00314876, -0.07473735, -0.99719828])
        # dCos = [-0.00314876, -0.07473735, -0.99719828]
        W = .55
        Telescopio.Trace(pSource_0, dCos, W)
        Rayos.push()







# ______________________________________#
#

Kos.display3d(Telescopio, Rayos ,1)

# Kos.display2d(Telescopio, Rayos,1,0)

# print(Telescopio.EFFL)
# X, Y, Z, L, M, N = Rayos.pick(-1)

# # ______________________________________#

# plt.plot(X, Y, 'x')
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Spot Diagram')
# plt.axis('square')
# plt.show()

