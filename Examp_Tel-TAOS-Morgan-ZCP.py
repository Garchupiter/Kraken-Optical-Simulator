# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

import KrakenOS as Kn
import numpy as np
import matplotlib.pyplot as plt

##############################################################
P_Obj = Kn.surf()

P_Obj.Thickness = 1772.336277098806
P_Obj.Glass = "AIR"
P_Obj.Diameter = 6.796727741707513E+002 * 2.0
P_Obj.Drawing = 0

M1 = Kn.surf()
M1.Rc = -6.06044E+003
M1.Thickness = (-1.774190000000000E+003) + (1.853722901194000E+000)
M1.k = -1.637E+000
M1.Glass = "MIRROR"
M1.Diameter = 6.63448E+002 * 2.0
M1.InDiameter = 228.6 * 2.0

#############################

M2 = Kn.surf()
M2.Rc = -6.06044E+003
M2.Thickness = -M1.Thickness
M2.k = -3.5782E+001
M2.Glass = "MIRROR"
M2.Diameter = 2.995730651164167E+002 * 2.0
ED0 = np.zeros(20)
ED0[2] = 4.458178314555000E-018
M2.AspherData = ED0
M2.AxisMove = 0

Vertex = Kn.surf()
Vertex.Thickness = 30.0
Vertex.Glass = "AIR"
Vertex.Diameter = 600.0
Vertex.Drawing = 0

Corrector_c1 = Kn.surf()
Corrector_c1.Thickness = 6.6E+000
Corrector_c1.Glass = "SILICASCHOTT"  # "BK7"#"LITHOSIL-Q"
Corrector_c1.Diameter = 118.0 * 2
ED1 = np.zeros(20)
ED1[0] = 2.059727552003000E-005
ED1[1] = -1.135080732384000E-009
Corrector_c1.AspherData = ED1

Corrector_c2 = Kn.surf()
Corrector_c2.Thickness = 341.65484183207997 - 100.0
Corrector_c2.Glass = "AIR"
Corrector_c2.Diameter = 118.0 * 2.0

Corrector_c2 = Kn.surf()
Corrector_c2.Thickness = 341.65484183207997
Corrector_c2.Glass = "AIR"
Corrector_c2.Diameter = 118.0 * 2

P_Ima = Kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 300.0
P_Ima.Drawing = 0
P_Ima.Name = "Plano imagen"

A = [P_Obj, M1, M2, Vertex, Corrector_c1, Corrector_c2, P_Ima]
configuracion_1 = Kn.Kraken_setup()

Telescopio = Kn.system(A, configuracion_1)

Rayos0 = Kn.raykeeper(Telescopio)

#######################################


W = 0.4
sup = 1
AperVal = 1300
AperType = "EPD"  # "STOP"
Pup = Kn.PupilCalc(Telescopio, sup, W, AperType, AperVal)
Pup.Samp = 5
Pup.Ptype = "hexapolar"
Pup.FieldType = "angle"
Pup.FieldX = 0.

Telescopio.IgnoreVignetting(0)

rays = Pup.Pattern2Field()

x, y, z, l, m, n = rays
W = 0.4
for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [l[i], m[i], n[i]]
    Telescopio.Trace(pSource_0, dCos, W)

    Rayos0.push()

############################################

Kn.display3d(Telescopio, Rayos0, 0)

X, Y, Z, L, M, N = Rayos0.pick(-1)

cenX = np.mean(X)
cenY = np.mean(Y)

plt.plot(X, Y, 'x')

# axis labeling
plt.xlabel('numbers')
plt.ylabel('values')

# figure name
plt.title('Dot Plot : Red Dots')
plt.axis('square')
plt.show()
