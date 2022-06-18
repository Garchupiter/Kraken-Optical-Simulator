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
        self.system.SDT[1].Rc = V[0]
        self.system.SDT[1].Thickness = V[1]
        self.system.SDT[2].Thickness = -V[1]
        self.system.SDT[2].Rc = V[0]
        self.system.SetData()

        self.system.Parax(self.w)
        D_EFFL = self.system.EFFL - self.effl

        Kos.TraceLoop(self.x, self.y, self.z, self.L, self.M, self.N, self.w, \
                      self.raykeeper, clean = 1)

        B = BestRMS_2(self.system, self.raykeeper)
        self.system.RestoreData()
        return [D_EFFL, B]


# ______________________________________#

currentDirectory = os.getcwd()
sys.path.insert(1, currentDirectory + '/library')

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Thickness = 2000.0
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1300.0
P_Obj.Drawing = 0

Thickness = 1778.0
M1 = Kos.surf()
M1.Rc = -6081.3
M1.Thickness = -Thickness
M1.k = -1.0
M1.Glass = "MIRROR"
M1.Diameter = 1300
M1.InDiameter = 650.0 * 0.5

M2 = Kos.surf()
M2.Rc = -6081.3
M2.Thickness = Thickness
M2.k = 0.0
M2.Glass = "MIRROR"
M2.Diameter = 650.0

M1Vertex=Kos.surf(Thickness = 50.0, Glass = "AIR",Diameter = 600.0)
M1Vertex.Drawing=0

Corr_1=Kos.surf(Thickness = 10, Glass = "LITHOSIL-Q", Diameter=230)

ED1=np.zeros(20)
ED1[0] = 0.0
ED1[1] = 0.0
Corr_1.AspherData = ED1

Corr_2=Kos.surf(Thickness=390, Glass="AIR", Diameter=230)
Corr_2.AspherData = -ED1

P_Ima = Kos.surf()
P_Ima.Diameter = 300.0

A = [P_Obj, M1, M2, M1Vertex, Corr_1, Corr_2, P_Ima]
configuracion_1 = Kos.Setup()
Telescope = Kos.system(A, configuracion_1)
# ______________________________________#


Rays = Kos.raykeeper(Telescope)


W = 0.6
effl = 5200
Surf=1

MyFun = Function2Optimize(Telescope, Rays, effl, Surf, W, M1.Diameter, "EPD")


R0 = -6081.3
D  = -1778.0


LimInf = [-6081.3 - 100, -1778.0 - 100]
LimSup = [-6081.3 + 100, -1778.0 + 100]
b=(LimInf, LimSup)

R = scipy.optimize.least_squares(MyFun.EFFL, [R0, D],bounds=b, verbose=2)

[R0, D] = R.x

Telescope.SDT[1].Rc = R0
Telescope.SDT[1].Thickness = D
Telescope.SDT[2].Thickness = -D
Telescope.SDT[2].Rc = R0
Telescope.SetData()
Telescope.SetSolid()


Surf, W, AperVal, AperType = 1, W, M1.Diameter, "EPD"
P = Kos.PupilCalc(Telescope, Surf, W, AperType, AperVal)
P.Samp, P.Ptype, P.FieldX, P.FieldType = 7, "fany", 0, "angle"
x, y, z, L, M, N = P.Pattern2Field()

Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)
x, y, z, l, m, n = Rays.pick(-1 , coordinates="local")

rms = R_RMS_delta(0,l,m,n,x,y)
print("R0: ", R0, " R1: ", -R0, " D: ", D, " RMS: ", rms )

Kos.display3d(Telescope, Rays, 2)



Telescope.Parax(W)
print(Telescope.EFFL)


P.Samp, P.Ptype, P.FieldX, P.FieldType = 7, "hexapolar", 0, "angle"
x, y, z, L, M, N = P.Pattern2Field()
Kos.TraceLoop(x, y, z, L, M, N, W, Rays, clean = 1)
RAYS = [Rays]
MK = ["x"]
COL = [[0.8,0.0,0.0]]
SURF = [-1]
MyPlot(RAYS, SURF, figure= "Spot2", mk = MK, col = COL)

