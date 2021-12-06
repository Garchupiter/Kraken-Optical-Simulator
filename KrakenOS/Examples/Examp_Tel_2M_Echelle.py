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
P_Obj.Thickness = 1000*0 + 3452.2
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1059. * 2.0
P_Obj.Drawing=0

# ______________________________________#

Thickness = 3452.2
M1 = Kos.surf()
M1.Rc = -9638.0
M1.Thickness = -Thickness
M1.k = -1.07731
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
M1.InDiameter = 250 * 2.0
M1.TiltY = 0.0
M1.TiltX = 0.0
M1.AxisMove = 0
# ______________________________________#

M2 = Kos.surf()
M2.Rc = -3930.0
M2.Thickness = Thickness + 1037.525880
M2.k = -4.3281
M2.Glass = "MIRROR"
M2.Diameter = 336.5 * 2.0
M2.TiltY = 0.0
M2.TiltX = 0.0
M2.DespY = 0.0
""" Se inclina el secundario """
M2.DespX = 1.0
M2.AxisMove = 0

# ______________________________________#


N1 = Kos.surf()
N1.Thickness = 720
N1.Glass = "NULL"
N1.Diameter = 100
N1.TiltX = 10
N1.Order = 0
N1.AxisMove = 1

############################

Colim = Kos.surf()
Colim.Rc = -1440.0
Colim.k = -1.0
Colim.Thickness = -720
Colim.Glass = "MIRROR"
Colim.Diameter = 400
Colim.TiltX = 0
Colim.Order = 0

############################
N2 = Kos.surf()
N2.Glass = "NULL"
N2.DespY =125.9836755
N2.AxisMove = 1
############################

N3 = Kos.surf()
N3.Glass = "NULL"
N3.TiltX =-71
N3.AxisMove = 1

#############################

G1 = Kos.surf()
G1.Glass = "MIRROR"
G1.Diameter = 300
G1.Grating_D = (1/0.079)
G1.Diff_Ord = -41
G1.Grating_Angle = 10.

#############################

N4 = Kos.surf()
N4.Glass = "NULL"
N4.TiltX = -N3.TiltX
N4.AxisMove = 1

#############################

N5 = Kos.surf()
N5.Glass = "NULL"
N5.Thickness = 1000
N5.TiltX = -N1.TiltX
N5.AxisMove = 1

#############################

# N6 = Kos.surf()
# N6.Glass = "NULL"
# N6.TiltX = -4.95736765
# N6.DespY = 47.70641809
# N6.AxisMove = 1

# #############################

# N7 = Kos.surf()
# N7.Glass = "NULL"
# N7.TiltZ = 90

# #############################

# N8 = Kos.surf()
# N8.Glass = "NULL"
# N8.TiltX = 20.7

# #############################

# G2 = Kos.surf()
# G2.Glass = "AIR"
# G2.Diameter = 18300
# G2.Thickness = 100
# # G2.Grating_D = 0#1/0.3
# # G2.Diff_Ord = 0
# G2.Grating_Angle = 0.0



P_Ima = Kos.surf()
P_Ima.Diameter = 1900
P_Ima.Rc = -1000.
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, M1, M2, N1, Colim, N2, N3, G1,N4, N5, P_Ima]
configuracion_1 = Kos.Setup()
Telescopio = Kos.system(A, configuracion_1)

# ______________________________________#
Telescopio = Kos.system(A, configuracion_1)
Rayos = Kos.raykeeper(Telescopio)

# ______________________________________#

tam = 9
rad = 2100 / 2
tsis = len(A) - 1

# ______________________________________#



# ______________________________________#

# for i in range(-tam, tam + 1):
for j in range(-tam, tam + 1):
    x_0 = (0 / tam) * rad
    y_0 = (j / tam) * rad
    r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
    if r < rad:
        print(0)
        tet = 0.0
        pSource_0 = [x_0, y_0, 0.0]
        dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
        W = .55
        Telescopio.Trace(pSource_0, dCos, W)
        Rayos.push()

# ______________________________________#


Kos.display2d(Telescopio, Rayos,0,0)

# Kos.display2d(Telescopio, Rayos,1,0)

print(Telescopio.EFFL)
X, Y, Z, L, M, N = Rayos.pick(-1)

# ______________________________________#

plt.plot(X, Y, 'x')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Spot Diagram')
plt.axis('square')
plt.show()

