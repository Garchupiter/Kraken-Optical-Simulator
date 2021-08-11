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
P_Obj.Thickness = 1000.0
P_Obj.Diameter = 300
P_Obj.Drawing = 0

M1 = kn.surf()
M1.Rc = -2000.0
M1.Thickness = M1.Rc / 2
M1.k = -1.0
M1.Glass = "MIRROR"
M1.Diameter = 300
M1.ShiftY = 200

P_Ima = kn.surf()
P_Ima.Glass = "AIR"
P_Ima.Diameter = 1600.0
P_Ima.Drawing = 0
P_Ima.Name = "Plano imagen"

A = [P_Obj, M1, P_Ima]

######################


configuracion_1 = kn.Kraken_setup()

Espejo = kn.system(A, configuracion_1)
Rayos = kn.raykeeper(Espejo)

tam = 5
rad = 150.0
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
            Espejo.Trace(pSource_0, dCos, W)
            Rayos.push()

            # W=0.5
            # Espejo.Trace(pSource_0, dCos,W)
            # Rayos.push()

            # W=0.6
            # Espejo.Trace(pSource_0, dCos,W)
            # Rayos.push()

kn.display3d(Espejo, Rayos, 0)
