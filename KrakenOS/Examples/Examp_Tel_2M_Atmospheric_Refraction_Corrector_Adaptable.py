# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp-Tel_2M_Atmospheric_Refraction_Corrector"""

import pkg_resources
required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("No instalado")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos
import numpy as np
import matplotlib.pyplot as plt
import scipy





def MyPlot(RK,surf, figure= "Spot", mk=["x"], col = [[0.8,0.0,0.0]]):
    fig = plt.figure(figure)
    ax = fig.add_subplot()
    i=0
    for rk in RK:
        x,y,z,l,m,n = rk.pick(surf[i], coordinates="local")
        plt.scatter(x, y,marker = mk[i], color = col[i])
        i = i + 1
    # https://matplotlib.org/3.1.1/api/markers_api.html
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('Spot diagram')
    plt.show()



def R_RMS_delta(Z1, L, M, N, X0, Y0):
    X1 = ((L / N) * Z1) + X0
    Y1 = ((M / N) * Z1) + Y0
    cenX = np.mean(X1)
    cenY = np.mean(Y1)
    x1 = (X1 - cenX)
    y1 = (Y1 - cenY)
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    return R_RMS

def BestFocus(X, Y, Z, L, M, N, system, mod=1):
    delta_Z = 0
    ZZ = (L, M, N, X, Y)
    v = scipy.optimize.fsolve(R_RMS_delta, delta_Z, args=ZZ)
    if mod ==1:
        system.SDT[-2].Thickness = system.SDT[-2].Thickness + v[0]
        system.SetData()
        system.SetSolid()
    return system, v[0]


def BestRMS(system,raykeeper):
    x,y,z,l,m,n = raykeeper.pick(-1 , coordinates="local")
    system, deltaZ = BestFocus(x, y, z, l, m, n, system, mod = 0)
    rms = R_RMS_delta(deltaZ,l,m,n,x,y)
    return rms
# _________________________________________________________________#

T = 0 # Prism rotation
ZenitDist = 55.0


A = 1.55/2 *0
P_Obj = Kos.surf(Thickness = 4452.2, Glass = "AIR", Diameter = 2118.0)
P_Obj.Drawing = 0

Thickness = 3452.2
M1 = Kos.surf()
M1.Rc = -9638.0
M1.Thickness = -Thickness
M1.k = -1.07731
M1.Glass = "MIRROR"
M1.Diameter = 1059.0 * 2.0
M1.InDiameter = 250 * 2.0

M2 = Kos.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness
M2.k = -4.3281
M2.Glass = "MIRROR"
M2.Diameter = 336.5 * 2.0

M1Vertex = Kos.surf()
M1Vertex.Thickness = 737.525
M1Vertex.Diameter = 200
M1Vertex.TiltZ=T
M1Vertex.AxisMove = 1
M1Vertex.Drawing =0


C1 = Kos.surf()
C1.Thickness = 5
C1.Glass = "BK7"
C1.Diameter = 100

C2 = Kos.surf()
C2.Thickness = C1.Thickness
C2.Glass = "F2"
C2.Diameter = 100
C2.TiltY = A
C2.AxisMove = 0

C3 = Kos.surf()
C3.Thickness = C1.Thickness
C3.Glass = "AIR"
C3.Diameter = 100
C3.TiltZ = -2 * M1Vertex.TiltZ
C3.AxisMove = 1


C4 = Kos.surf()
C4.Thickness = C1.Thickness
C4.Glass = "BK7"
C4.Diameter = 100

C5 = Kos.surf()
C5.Thickness = C1.Thickness
C5.Glass = "F2"
C5.Diameter = 100
C5.TiltY = A
C5.AxisMove = 0

C6 = Kos.surf()
C6.Thickness = C1.Thickness + 277.252
C6.Glass = "AIR"
C6.Diameter = 100

RZ = Kos.surf(Diameter = 100)
RZ.TiltZ =  M1Vertex.TiltZ
RZ.AxisMove = 1

P_Ima = Kos.surf(Diameter = 100.0)

A = [P_Obj, M1, M2, M1Vertex, C1, C2, C3, C4, C5, C6, RZ, P_Ima]
Config_1 = Kos.Setup()
Tel = Kos.system(A, Config_1)

# _________________________________________________________________#

Rays1 = Kos.raykeeper(Tel)
Rays2 = Kos.raykeeper(Tel)
Rays3 = Kos.raykeeper(Tel)





W = 0.60169
sup = 1  # Difining M1 as enter pupil diameter
AperVal = 2000
AperType = "EPD"  # "STOP"
Pup = Kos.PupilCalc(Tel, sup, W, AperType, AperVal)
Pup.Samp = 11
Pup.FieldType = "angle"

Pup.AtmosRef = 1
Pup.T = 283.15  # k
Pup.P = 101300  # Pa
Pup.H = 0.5  # Humidity ratio 1 to 0
Pup.xc = 400  # ppm
Pup.lat = 31  # degrees
Pup.h = 2800  # meters
Pup.l1 = W  # micron
Pup.l2 = 0.50169  # micron
Pup.z0 = ZenitDist  # degrees
Pup.Ptype = "hexapolar"
Pup.FieldX = 0.0


W = [0.50169, 0.60169, 0.70169]
RAYS=[Rays1, Rays2, Rays3]
i = 0
for w in W:
    Pup.l2 = w
    x, y, z, L, M, N = Pup.Pattern2Field()
    Kos.TraceLoop(x, y, z, L, M, N, w, RAYS[i], clean = 1)
    i = i + 1

MK = [".", ".", "."]
COL = [[0.8,0.0,0.0], [0.0,0.8,0.0], [0.0,0.0,0.8]]
SURF = [-1, -1, -1]
MyPlot(RAYS, SURF, figure= "Spot", mk = MK, col = COL)
