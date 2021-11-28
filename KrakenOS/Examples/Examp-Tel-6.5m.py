
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

from scipy import optimize
from scipy.optimize import fmin_cg

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

import numpy as np
import matplotlib.pyplot as plt
import time


P_Obj=Kos.surf()
P_Obj.Rc=0
P_Obj.Thickness=3000+3.452200000000000E+003
P_Obj.Glass="AIR"
P_Obj.Diameter=6502.4
P_Obj.Drawing=0

Thickness=6178
M1=Kos.surf()
M1.Rc=-16256.0
M1.Thickness=-Thickness
M1.k=-1.
M1.Glass="MIRROR"
M1.Diameter=6502.4

M2=Kos.surf()
M2.Rc=-5150.974
M2.Thickness=Thickness+1851.28
M2.k=-2.6946
M2.Glass="MIRROR"
M2.Diameter=1714.5
M2.DespX=0
M2.AxisMove=0

P_Ima=Kos.surf()
P_Ima.Diameter=1000.0
P_Ima.Glass="AIR"
P_Ima.Name="Plano imagen"

A = [P_Obj, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()
Telescope = Kos.system(A, configuracion_1)

# ______________________________________#

W = 0.4861
sup = 1
AperVal = 6500
AperType = "STOP"  # "STOP"
Pup = Kos.PupilCalc(Telescope, sup, W, AperType, AperVal)

Rays = Kos.raykeeper(Telescope)

# _________________________________________#

print("Radio pupila de entrada: ")
print(Pup.RadPupInp)
print("Posición pupila de entrada: ")
print(Pup.PosPupInp)
print("Rádio pupila de salida: ")
print(Pup.RadPupOut)
print("Posicion pupila de salida: ")
print(Pup.PosPupOut)
print("Posicion pupila de salida respecto al plano focal: ")
print(Pup.PosPupOutFoc)
print("Orientación pupila de salida")
print(Pup.DirPupSal)

print("Airy disk diameter focal distance (micrometers)")
print(Pup.FocusAiryRadius)


print("EFFL: ", Telescope.EFFL)
# _________________________________________#

[L, M, N] = Pup.DirPupSal
print(L, M, N)




Pup.Samp = 1
Pup.Ptype = "chief"
Pup.FieldType = "angle"
Pup.FieldX = 0.0
x, y, z, L, M, N = Pup.Pattern2Field()

pSource_0 = [x[0]+6500/2, y[0], z[0]]
dCos = [L[0], M[0], N[0]]
Telescope.Trace(pSource_0, dCos, W)
Rays.push()



print(" ")
print(" Coordenadas ")
print(Telescope.XYZ[0])
print(Telescope.XYZ[1])
print(Telescope.XYZ[2])
print(Telescope.XYZ[3][0])



print(" ")
print(" Cosenos directores ")
print(Telescope.LMN[0])
print(Telescope.LMN[1])
print(Telescope.LMN[2])





print(" ")
print(" Distancias ")
print(Telescope.DISTANCE[0])
print(Telescope.DISTANCE[1])
print(Telescope.DISTANCE[2])



# print(" ")
# print(" op ")
# print(" ")
# print(" Distancias opticas")
# print(Telescope.OP[0])
# print(Telescope.OP[1])
# print(Telescope.OP[2])
# print(Telescope.OP[3])
# print(Telescope.OP[4])
# print(Telescope.OP[5])
# print(Telescope.OP[6])
# print(Telescope.OP[7])
# print(Telescope.OP[8])
# print(Telescope.OP[9])
# print("Total = ", np.sum(Telescope.OP)- np.sum(Telescope.OP[0]))
# print("    ")

# print(Telescope.RP)
# print(Telescope.RS)
# print(Telescope.TP)
# print(Telescope.TS)
# print("Total trans : ")
# print(Telescope.TT)


















W = 0.4861
sup = 1
AperVal = 6500.
AperType = "EPD"
Pup = Kos.PupilCalc(Telescope, sup, W, AperType, AperVal)

# _________________________________________#

print("Radio pupila de entrada: ")
print(Pup.RadPupInp)
print("Posición pupila de entrada: ")
print(Pup.PosPupInp)
print("Rádio pupila de salida: ")
print(Pup.RadPupOut)
print("Posicion pupila de salida: ")
print(Pup.PosPupOut)
print("Posicion pupila de salida respecto al plano focal: ")
print(Pup.PosPupOutFoc)
print("Orientación pupila de salida")
print(Pup.DirPupSal)

print("Airy disk diameter focal distance (micrometers)")
print(Pup.FocusAiryRadius)


print("EFFL: ", Telescope.EFFL)
# _________________________________________#

[L, M, N] = Pup.DirPupSal
print(L, M, N)




# Pup.Samp = 1
# Pup.Ptype = "chief"
# Pup.FieldType = "angle"
# Pup.FieldX = 1.0
# x, y, z, L, M, N = Pup.Pattern2Field()

# pSource_0 = [x[0], y[0], z[0]]
# dCos = [L[0], M[0], N[0]]
# Telescope.Trace(pSource_0, dCos, W)
# Rays.push()



# print(" ")
# print(" Coordenadas ")
# print(Telescope.XYZ[0])
# print(Telescope.XYZ[1])
# print(Telescope.XYZ[2])
# print(Telescope.XYZ[3])
# print(Telescope.XYZ[4])
# print(Telescope.XYZ[5])
# print(Telescope.XYZ[6])
# print(Telescope.XYZ[7])
# print(Telescope.XYZ[8])
# print(Telescope.XYZ[9])
# XYZ= Telescope.XYZ[10]
# print(XYZ[0], XYZ[1], XYZ[2])


# print(" ")
# print(" Cosenos directores ")
# print(Telescope.LMN[0])
# print(Telescope.LMN[1])
# print(Telescope.LMN[2])
# print(Telescope.LMN[3])
# print(Telescope.LMN[4])
# print(Telescope.LMN[5])
# print(Telescope.LMN[6])
# print(Telescope.LMN[7])
# print(Telescope.LMN[8])
# print(Telescope.LMN[9])
# LMN = Telescope.LMN[9]
# print(LMN[0], LMN[1], LMN[2])





# print(" ")
# print(" Distancias ")
# print(Telescope.DISTANCE[0])
# print(Telescope.DISTANCE[1])
# print(Telescope.DISTANCE[2])
# print(Telescope.DISTANCE[3])
# print(Telescope.DISTANCE[4])
# print(Telescope.DISTANCE[5])
# print(Telescope.DISTANCE[6])
# print(Telescope.DISTANCE[7])
# print(Telescope.DISTANCE[8])
# print(Telescope.DISTANCE[9])

# print(" ")
# print(" op ")
# print(" ")
# print(" Distancias opticas")
# print(Telescope.OP[0])
# print(Telescope.OP[1])
# print(Telescope.OP[2])
# print(Telescope.OP[3])
# print(Telescope.OP[4])
# print(Telescope.OP[5])
# print(Telescope.OP[6])
# print(Telescope.OP[7])
# print(Telescope.OP[8])
# print(Telescope.OP[9])
# print("Total = ", np.sum(Telescope.OP)- np.sum(Telescope.OP[0]))
# print("    ")

# print(Telescope.RP)
# print(Telescope.RS)
# print(Telescope.TP)
# print(Telescope.TS)
# print("Total trans : ")
# print(Telescope.TT)


# _________________________________________#

Pup.Samp = 15
Pup.Ptype = "hexapolar"
Pup.FieldType = "angle"
Pup.FieldX =0.5

# x, y, z, L, M, N = Pup.Pattern2Field()
# for i in range(0, len(x)):
#     pSource_0 = [x[i], y[i], z[i]]
#     dCos = [L[i], M[i], N[i]]
#     Telescope.Trace(pSource_0, dCos, W)
#     Rays.push()

# Kos.display2d(Telescope, Rays, 1, 0 )



X, Y, Z, P2V = Kos.Phase(Pup)


"""Indicamos el grado de expanción para los polinomios de Zernike"""
NC = 35

"""Generamos un arreglo numpy conlas mismas dimensiones de la expanción definida"""
A = np.ones(NC)

"""Calculamos los polinomios de Zernike con la fase calculada y el numero de
  elementos deseados en la expanción, Zcoef son los coeficientes en longitudes de
  onda, Mat es la expreción matematica de Zeidel para dicho coeficiente,
  esto con fines ilustrativos, w_rms es el error del ajuste"""

Zcoef, Mat, RMS2Chief, RMS2Centroid, FITTINGERROR = Kos.Zernike_Fitting(X, Y, Z, A)

"""" Se despliegan los resultados """

for i in range(0, NC):
    print("z", i + 1, "  ", "{0:.8f}".format(float(Zcoef[i])), ":", Mat[i])

# # ______________________________________#
print("(RMS) Fitting error: ", FITTINGERROR)
print(RMS2Chief, "RMS(to chief) From fitted coefficents")
print(RMS2Centroid, "RMS(to centroid) From fitted coefficents")
