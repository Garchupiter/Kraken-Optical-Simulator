#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Perfect Lens"""

import time
import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("No instalado")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos
from scipy.optimize import fsolve
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


class Function2Optimize:
    def __init__(self, system, raykeeper, effl, Surf, w, AperVal, AperType):

        self.system = system
        self.w = w
        self.effl = effl
        self.raykeeper = raykeeper

        self.P = Kos.PupilCalc(Lens, Surf, w, AperType, AperVal)

        self.P.Samp, self.P.Ptype, self.P.FieldType = 5, "fanx", "angle"
        self.P.FieldX, self.P.FieldY = 0.0, 0.0
        self.ProcessPattern2Field()

    def ProcessPattern2Field(self):
        self.x, self.y, self.z, self.L, self.M, self.N = self.P.Pattern2Field()

    def EFFL(self, R):
        self.system.SDT[1].Rc = R[0]
        self.system.SDT[2].Rc = -R[0]
        self.system.SetData()
        self.system.Parax(self.w)
        Kos.TraceLoop(self.x, self.y, self.z, self.L, self.M, self.N, self.w, self.raykeeper, clean = 1)

        D_EFFL = self.system.EFFL - self.effl

        # B = BestRMS(self.system, self.raykeeper)

        self.system.RestoreData()
        return [D_EFFL]




P_Obj = Kos.surf()
P_Obj.Thickness = 100
P_Obj.Diameter = 120.0

Fa = Kos.surf()
Fa.Rc = 0.0
Fa.Thickness = 12.0
Fa.Glass = "BK7"
Fa.Diameter = 120.0

Fb = Kos.surf()
Fb.Rc = 0.0
Fb.Thickness = 1000.0
Fb.Diameter = 120
Fb.Glass = "AIR"


P_Ima = Kos.surf()
P_Ima.Thickness = 0.0
P_Ima.Diameter = 120.0
P_Ima.Name = "Image plane"
P_Ima.Nm_Pos = (-10,10)
P_Ima.Glass = "AIR"

config_1 = Kos.Setup()

A = [P_Obj, Fa, Fb, P_Ima]
Lens = Kos.system(A, config_1)

Rays = Kos.raykeeper(Lens)







Wave = 0.5
effl = 1000
Lens.Parax(0.5)
Surf=1

print("Initial effective focal length: ", Lens.EFFL)
# print("Initial RMS radius in best focus: ", BestRMS(Lens, Wave))

MyFun = Function2Optimize(Lens, Rays, effl, Surf, Wave, Fa.Diameter, "EPD")

R0 = -1000000
R1 = 1000000
R0 = fsolve(MyFun.EFFL, [R0])

print("Resulting R0: ", R0, " R1: ", R0)
Lens.SDT[1].Rc = R0
Lens.SDT[2].Rc = -R0
Lens.SetData()
Lens.SetSolid()


Surf, W, AperVal, AperType = 1, 0.45, Fa.Diameter, "EPD"
P = Kos.PupilCalc(Lens, Surf, W, AperType, AperVal)
P.Samp, P.Ptype, P.FieldX, P.FieldType = 8, "fanx", 0, "angle"
x, y, z, L, M, N = P.Pattern2Field()

Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)


Lens.Parax(0.5)
print("Final effective focal length: ", Lens.EFFL)



x, y, z, l, m, n = Rays.pick(-1 , coordinates="local")
system, deltaZ = BestFocus(x, y, z, l, m, n, Lens, mod = 1)

x, y, z, L, M, N = P.Pattern2Field()
Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)

x, y, z, l, m, n = Rays.pick(-1 , coordinates="local")
rms = R_RMS_delta(deltaZ,l,m,n,x,y)
print("Final RMS radius: ", rms)




Kos.display2d(Lens, Rays, 1)



P.Samp, P.Ptype, P.FieldX, P.FieldType = 7, "hexapolar", 0, "angle"
x, y, z, L, M, N = P.Pattern2Field()
Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)
RAYS = [Rays]
MK = ["x"]
COL = [[0.8,0.0,0.0]]
SURF = [-1]
MyPlot(RAYS, SURF, figure= "Spot2", mk = MK, col = COL)

