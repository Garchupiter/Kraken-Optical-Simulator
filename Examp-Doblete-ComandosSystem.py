#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

import numpy as np

import Kraken as kn

##############################################################    
P_Obj = kn.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 3.452229924716749E+003
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1.059E+003 * 2.0
P_Obj.Drawing = 0

M1 = kn.surf()
M1.Rc = -9.638000000004009E+003
M1.Thickness = -P_Obj.Thickness
M1.k = -1.077310000000000E+000
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
M1.InDiameter = 250 * 2.0
M1.TiltX = 0.0
M1.AxisMove = 0

M2 = kn.surf()
M2.Rc = -3.93E+003
M2.Thickness = P_Obj.Thickness + 1.037179115116706E+003
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0

P_Ima = kn.surf()
P_Ima.Diameter = 100.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

A = [P_Obj, M1, M2, P_Ima]

######################


configuracion_1 = kn.Kraken_setup()

Doblete = kn.system(A, configuracion_1)
Rayos = kn.raykeeper(Doblete)

tam = 10
rad = 1000.0
tsis = len(A) - 1

for i in range(-tam, tam + 1):
    for j in range(-tam, tam + 1):

        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]

            W = 0.4
            Doblete.Trace(pSource_0, dCos, W)
            Rayos.push()

print("Distancia focal efectiva")
print(Doblete.EFFL)
print("Plano principal anterior")
print(Doblete.PPA)
print("Plano principal posterior")
print(Doblete.PPP)
print("Superficies tocadas por el rayo")
print(Doblete.SURFACE)
print("Nombre de la superficie")
print(Doblete.NAME)
print("Vidrio de la superficie")
print(Doblete.GLASS)
print("Coordenadas del rayo en las superficies")
print(Doblete.XYZ)
print("Etc, ver documentaciòn")
print(Doblete.S_XYZ)
print(Doblete.T_XYZ)
print(Doblete.OST_XYZ)
print(Doblete.DISTANCE)
print(Doblete.OP)
print(Doblete.TOP)
print(Doblete.TOP_S)
print(Doblete.ALPHA)
print(Doblete.BULK_TRANS)
print(Doblete.S_LMN)
print(Doblete.LMN)
print(Doblete.R_LMN)
print(Doblete.N0)
print(Doblete.N1)
print(Doblete.WAV)
print(Doblete.G_LMN)
print(Doblete.ORDER)
print(Doblete.GRATING)
print(Doblete.RP)
print(Doblete.RS)
print(Doblete.TP)
print(Doblete.TS)
print(Doblete.TTBE)
print(Doblete.TT)
