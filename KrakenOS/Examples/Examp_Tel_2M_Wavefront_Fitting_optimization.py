# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M Wavefront Fitting"""

import os
import sys
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

# ______________________________________#

currentDirectory = os.getcwd()
sys.path.insert(1, currentDirectory + '/library')

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 1000 + 3.452200000000000E+003
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1.059E+003 * 2.0

# ______________________________________#

Thickness = 3.452200000000000E+003
M1 = Kos.surf()
M1.Rc = -9.638000000004009E+003
M1.Thickness = -Thickness
M1.k = -1.077310000000000E+000
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
M1.InDiameter = 250 * 2.0
M1.TiltY = 0.0
M1.TiltX = 0.0

# ______________________________________#

M1.AxisMove = 0
M2 = Kos.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness + 1037.525880
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0
M2.TiltY = 0.0
M2.TiltX = 0.0
M2.DespY = 0.0
M2.DespX = 0.0
M2.AxisMove = 0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 300.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()
Telescopio = Kos.system(A, configuracion_1)

# ______________________________________#

Surf = 1
W = 0.5016
AperVal = 2000.
AperType = "EPD"
Pupil = Kos.PupilCalc(Telescopio, Surf, W, AperType, AperVal)
Pupil.Samp = 10
Pupil.Ptype = "hexapolar"
Pupil.FieldX = 0.1
Pupil.FieldY = 0.0

Pupil.FieldType = "angle"

X, Y, Z, P2V = Kos.Phase(Pupil)
print("Peak to valley: ", P2V)
NC = 18
A = np.ones(38)

Zcoef, Mat, RMS2Chief, RMS2Centroid, FITTINGERROR = Kos.Zernike_Fitting(X, Y, Z, A)


for i in range(0, NC):
    print("z", i + 1, "  ", "{0:.6f}".format(float(Zcoef[i])), ":", Mat[i])
# ______________________________________#

# print("RMS: ", "{:.4f}".format(float(w_rms)), " Error del ajuste: ", fitt_error)
# z_coeff[0] = 0
# print("RMS to chief: ", np.sqrt(np.sum(z_coeff * z_coeff)))
# z_coeff[1] = 0
# z_coeff[2] = 0
# print("RMS to centroid: ", np.sqrt(np.sum(z_coeff * z_coeff)))

# #______________________________________#

RR = Kos.raykeeper(Telescopio)
x, y, z, L, M, N = Pupil.Pattern2Field()

# ______________________________________#

for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    RR.push()

# ______________________________________#

Kos.display2d(Telescopio, RR, 1,0)
X, Y, Z, L, M, N = RR.pick(-1)

# ______________________________________#

plt.plot(X, Y, 'x')
plt.xlabel('numbers')
plt.ylabel('values')
plt.title('spot Diagram')
plt.axis('square')
plt.show()

ima = Kos.WavefrontData2Image(Zcoef, 400)

Type = "interferogram"
Kos.ZernikeDataImage2Plot(ima, Type)



Surf = 1
W = 0.5016
AperVal = 2000.
AperType = "STOP"
Pupil = Kos.PupilCalc(Telescopio, Surf, W, AperType, AperVal)
Pupil.Samp = 10
Pupil.Ptype = "hexapolar"
Pupil.FieldX = 0.1
Pupil.FieldY = 0.0


AB = Kos.Seidel(Pupil)

print("--------------------------------------")
print(AB.SCW_AN)
print(AB.SCW_NM)
print(AB.SCW_TOTAL)


class Function2Optimize:
    def __init__(self, ABER):
        self.ABER = ABER

    def Fun(self, K):
        System = self.ABER.SYSTEM
        System.SDT[1].k = K[0]
        System.SDT[2].k = K[1]
        System.SetData()

        self.ABER.calculate()

        System.RestoreData()
        Sph = self.ABER.SAC_TOTAL[0]
        Coma = self.ABER.SAC_TOTAL[1]
        return [Sph, Coma]


SeidelFun = Function2Optimize(AB)

Pupil.FieldY = 0.2

from scipy.optimize import fsolve

root = fsolve(SeidelFun.Fun, [0, 0])
print(root)
