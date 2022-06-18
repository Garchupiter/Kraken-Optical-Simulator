# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M Wavefront Fitting"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
""" Looking for if KrakenOS is installed, if not, it assumes that
an library downloaded from github is running"""

required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("Not installed")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos

import scipy



###########-----------------------




import matplotlib.pyplot as plt
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


def BestRMS_2(system,raykeeper):
    x,y,z,l,m,n = raykeeper.pick(-1 , coordinates="local")
    rms = R_RMS_delta(0,l,m,n,x,y)
    return rms


def BestFocus(X, Y, Z, L, M, N, system, mod=1):
    delta_Z = 0
    ZZ = (L, M, N, X, Y)
    v = scipy.optimize.fsolve(R_RMS_delta, delta_Z, args=ZZ)
    if mod ==1:
        system.SDT[-2].Thickness = system.SDT[-2].Thickness + v[0]
        system.SetData()
        system.SetSolid()
    return system, v[0]

class Function2Optimize:
    def __init__(self, system, raykeeper, effl, Surf, w, AperVal, AperType):

        self.system = system
        self.w = w
        self.effl = effl
        self.raykeeper = raykeeper

        self.P = Kos.PupilCalc(system, Surf, w, AperType, AperVal)

        self.P.Samp, self.P.Ptype, self.P.FieldType = 25, "fanx", "angle"
        self.P.FieldX, self.P.FieldY = 0.0, 0.0
        self.ProcessPattern2Field()

    def ProcessPattern2Field(self):
        self.x, self.y, self.z, self.L, self.M, self.N = self.P.Pattern2Field()

    def EFFL(self, V):
        self.system.SDT[1].k = V[0]
        self.system.SDT[2].k = V[1]

        self.system.Parax(self.w)
        D_EFFL = self.system.EFFL - self.effl

        Kos.TraceLoop(self.x, self.y, self.z, self.L, self.M, self.N, self.w, \
                      self.raykeeper, clean = 1)

        B = BestRMS_2(self.system, self.raykeeper)
        self.system.RestoreData()
        return [D_EFFL, B]





class Function2OptimizeSeidel:
    def __init__(self, ABER, w, effl):
        self.ABER = ABER
        self.w =w
        self.effl = effl

    def Fun(self, V):
        System = self.ABER.SYSTEM
        System.SDT[1].k = V[0]
        System.SDT[2].k = V[1]
        ED1=np.zeros(20)
        ED1[0] = V[2]
        ED1[1] = V[3]

        System.SDT[4].AspherData =  ED1
        System.SDT[5].AspherData = -ED1



        System.SDT[1].Rc = V[4]
        System.SDT[2].Rc = V[4]

        System.SetData()

        System.Parax(self.w)
        D_EFFL = np.abs(System.EFFL - self.effl)

        self.ABER.calculate()

        System.RestoreData()
        Sph = self.ABER.SAC_TOTAL[0]
        Coma = self.ABER.SAC_TOTAL[1]
        Ast = self.ABER.SAC_TOTAL[2]
        Petz = self.ABER.SAC_TOTAL[3]

        return [Sph, Coma, Ast, Petz, D_EFFL]
# ______________________________________#

currentDirectory = os.getcwd()
sys.path.insert(1, currentDirectory + '/library')

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Thickness = 2000.0
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1300.0
P_Obj.Drawing = 0

Thickness = 1750.6548739435177
M1 = Kos.surf()
M1.Rc = -6034.3700365797185
M1.Thickness = -Thickness
M1.k = -1.0
M1.Glass = "MIRROR"
M1.Diameter = 1300
M1.InDiameter = 650.0 * 0.5

M2 = Kos.surf()
M2.Rc = -6034.3700365797185
M2.Thickness = Thickness
M2.k = 0.0
M2.Glass = "MIRROR"
M2.Diameter = 650.0

M1Vertex=Kos.surf(Thickness = 20.0, Glass = "AIR",Diameter = 600.0)
M1Vertex.Drawing=0

Corr_1=Kos.surf(Thickness = 10, Glass = "LITHOSIL-Q", Diameter=230)

ED1=np.zeros(20)
ED1[0] = 0.0
ED1[1] = 0.0
Corr_1.AspherData = ED1

Corr_2=Kos.surf(Thickness = 420, Glass="AIR", Diameter=230)
Corr_2.AspherData = -ED1

P_Ima = Kos.surf()
P_Ima.Diameter = 700.0

A = [P_Obj, M1, M2, M1Vertex, Corr_1, Corr_2, P_Ima]
configuracion_1 = Kos.Setup()
Telescope = Kos.system(A, configuracion_1)
Telescope.IgnoreVignetting(Prep=0)

# ______________________________________#


Rays = Kos.raykeeper(Telescope)
Rays0 = Kos.raykeeper(Telescope)
Rays1 = Kos.raykeeper(Telescope)
Rays2 = Kos.raykeeper(Telescope)


W = 0.6
Surf=1
AperVal = 1300.
AperType = "STOP"
Pupil = Kos.PupilCalc(Telescope, Surf, W, AperType, AperVal)
Pupil.Samp = 10
Pupil.Ptype = "hexapolar"
Pupil.FieldX = 1.7/2.0
Pupil.FieldY = 0.0

AB = Kos.Seidel(Pupil)
print("--------------------------------------")
print(AB.SCW_AN)
print(AB.SCW_NM)
print(AB.SCW_TOTAL)

effl = 5200

MyFun = Function2OptimizeSeidel(AB, W, effl)

K1 = -1.4855057854211016
K2 = -29.825480470600862
E1 = 0
E2 = 0
R1 =-6034.3700365797185


LimInf = [K1 - 3, K2 - 10, E1 - 1.e-8, E2 - 1.e-9, R1 - 100]
LimSup = [K1 + 3, K2 + 10, E1 + 1.e-8, E2 + 1.e-9, R1 + 100]
b=(LimInf, LimSup)

R = scipy.optimize.least_squares(MyFun.Fun, [K1, K2, E1, E2, R1],bounds=b, verbose=2)

[K1, K2, E1, E2, R1] = R.x

Telescope.SDT[1].k = K1
Telescope.SDT[2].k = K2

ED1=np.zeros(20)
ED1[0] = E1
ED1[1] = E2
Telescope.SDT[4].AspherData =  ED1
Telescope.SDT[5].AspherData = -ED1

Telescope.SDT[1].Rc = R1
Telescope.SDT[2].Rc = R1


Telescope.SetData()
Telescope.SetSolid()


Surf, W, AperVal, AperType = 1, W, M1.Diameter, "EPD"
P = Kos.PupilCalc(Telescope, Surf, W, AperType, AperVal)
P.Samp, P.Ptype, P.FieldX, P.FieldType = 11, "hexapolar", (1.7 / 2.0) , "angle"
x, y, z, L, M, N = P.Pattern2Field()

DeltaW = 0.1
Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)
Kos.TraceLoop(x, y, z, L, M, N, W + DeltaW, Rays, clean = 0)
Kos.TraceLoop(x, y, z, L, M, N, W - DeltaW, Rays, clean = 0)
x1, y1, z1, l1, m1, n1 = Rays.pick(-1 , coordinates="local")

BestFocus(x1,y1,z1,l1,m1,n1, Telescope, mod=1)





""" -----------------------------------------"""



P.Samp, P.Ptype, P.FieldX, P.FieldType = 7, "hexapolar", 0, "angle"
x, y, z, L, M, N = P.Pattern2Field()


Kos.TraceLoop(x, y, z, L, M, N, W, Rays0, clean = 1)
Kos.TraceLoop(x, y, z, L, M, N, W + DeltaW, Rays1, clean = 0)
Kos.TraceLoop(x, y, z, L, M, N, W - DeltaW, Rays2, clean = 0)
x, y, z, l, m, n = Rays0.pick(-1 , coordinates="local")


rms = R_RMS_delta(0,l,m,n,x,y)
print("K1: ", K1, " K2: ", K2, " E1: ", E1, " E2: ", E2, " RMS: ", rms )

# Kos.display3d(Telescope, Rays, 2)


RAYS = [Rays0, Rays1, Rays2]
MK = ["x", "*", "+"]
COL = [[0.8,0.0,0.0], [0.0,0.8,0.0], [0.0,0.0,0.8]]
SURF = [-1, -1, -1]
MyPlot(RAYS, SURF, figure= "Spot1", mk = MK, col = COL)


""" ---------------------------------------------"""

P.Samp, P.Ptype, P.FieldX, P.FieldType = 7, "hexapolar", 1.7/2, "angle"
x, y, z, L, M, N = P.Pattern2Field()

Kos.TraceLoop(x, y, z, L, M, N, W, Rays0, clean = 1)
Kos.TraceLoop(x, y, z, L, M, N, W + DeltaW, Rays1, clean = 1)
Kos.TraceLoop(x, y, z, L, M, N, W - DeltaW, Rays2, clean = 1)
x, y, z, l, m, n = Rays0.pick(-1 , coordinates="local")
system, deltaZ = BestFocus(x, y, z, l, m, n, Telescope, mod = 1)

rms = R_RMS_delta(0,l,m,n,x,y)
print("K1: ", K1, " K2: ", K2, " E1: ", E1, " E2: ", E2, " RMS: ", rms )


RAYS = [Rays0, Rays1, Rays2]
MK = ["x", "*", "+"]
COL = [[0.8,0.0,0.0], [0.0,0.8,0.0], [0.0,0.0,0.8]]
SURF = [-1, -1, -1]
MyPlot(RAYS, SURF, figure= "Spot2", mk = MK, col = COL)


Telescope.Parax(W)
print(Telescope.EFFL)

Kos.display2d(Telescope, Rays, 1)
Telescope.Parax(W)
print(Telescope.EFFL)













W = 0.6
Surf=1
AperVal = 1300.
AperType = "STOP"
Pupil = Kos.PupilCalc(Telescope, Surf, W, AperType, AperVal)
Pupil.Samp = 10
Pupil.Ptype = "hexapolar"
Pupil.FieldX = 1.7/2.0
Pupil.FieldY = 0.0

AB = Kos.Seidel(Pupil)
print("--------------------------------------")
print(AB.SCW_AN)
print(AB.SCW_NM)
print(AB.SCW_TOTAL)
