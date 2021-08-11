#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""
import Kraken as kn
import numpy as np

#############################################################


##############################################################    
P_Obj = kn.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 100
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0
P_Obj.Name = "P Obj"

L1a = kn.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "N-BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

L1b = kn.surf()
L1b.Rc = -3.071608670000159E+001

L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

L1c = kn.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001 - 40
L1c.Glass = "AIR"
L1c.Diameter = 30

pupila = kn.surf()
pupila.Rc = 0
pupila.Thickness = 40.
pupila.Glass = "AIR"
pupila.Diameter = 15.0
pupila.Name = "Ap Stop"

P_Ima = kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 20.0
P_Ima.Name = "Plano imagen"

A = [P_Obj, L1a, L1b, L1c, pupila, P_Ima]

##########################################################################


config_1 = kn.Kraken_setup()
Doblete = kn.system(A, config_1)

W = 0.6
sup = 4
AperVal = 3
AperType = "EPD"
field = 3.25
fieldType = "angle"

AB = kn.Seidel(Doblete, sup, W, AperType, AperVal, field, fieldType)
print( AB[0][0])
print(np.sum(AB[1][0]), np.sum(AB[1][1]), np.sum(AB[1][2]), np.sum(AB[1][3]), np.sum(AB[1][4]))

j=1
print( AB[0][0+j])
print(np.sum(AB[1+j][0]), np.sum(AB[1+j][1]), np.sum(AB[1+j][2]), np.sum(AB[1+j][3]), np.sum(AB[1+j][4]))

j=2
print( AB[0][0+j])
print(np.sum(AB[1+j][0]), np.sum(AB[1+j][1]), np.sum(AB[1+j][2]), np.sum(AB[1+j][3]), np.sum(AB[1+j][4]))

j=3
print( AB[0][0+j])
print(np.sum(AB[1+j][0]), np.sum(AB[1+j][1]), np.sum(AB[1+j][2]), np.sum(AB[1+j][3]), np.sum(AB[1+j][4]))






"""  kn.Seidel OUTPUT

Seidel Aberration Coefficents
SAC=[si, sii, siii, siv, sv]

Seidel coefficients in waves
SCW=[W040, W131, W222, W220, W311]

    
Transverse Aberration Coefficents
Spherical
TSPH = si / (2.0 * u[-1] * n_1[-1])
Coma
TSCO = sii / (2.0 * u[-1] * n_1[-1])
TTCO = 3.0 * sii / (2.0 * u[-1] * n_1[-1])
Astigmatism
TAST = -siii / (u[-1] * n_1[-1])
Field Curvature
TPFC = -siv / (2.0 * u[-1] * n_1[-1])
TSFC = (siii + siv) / (2.0 * u[-1] * n_1[-1])
TTFC = ((3.0 * siii) + siv) / (2.0 * u[-1] * n_1[-1])
Distortion
TDIS = -sv / (2.0 * u[-1] * n_1[-1])
TAC=[TSPH, TSCO, TTCO, TAST, TPFC, TSFC, TTFC, TDIS]
    

Longitudinal Aberration Coefficents
Spherical
LSPH = si / (2.0 * u[-1] * u[-1] * n_1[-1])
Coma
LSCO = sii / (2.0 * u[-1] * u[-1] * n_1[-1])
LTCO = 3.0 * sii / (2.0 * u[-1] * u[-1] * n_1[-1])
Astigmatism
LAST = siii / (u[-1] * u[-1] * n_1[-1])
Field Curvature
LPFC = siv / (2.0 * u[-1] * u[-1] * n_1[-1])
LSFC = (siii + siv) / (2.0 * u[-1] * u[-1] * n_1[-1])
LTFC = ((3.0 * siii) + siv) / (2.0 * u[-1] * u[-1] * n_1[-1])
## Distortion
LDIS = -sv / (2.0 * u[-1] * u[-1] * n_1[-1])
LAC=[LSPH, LSCO, LTCO, LAST, LPFC, LSFC, LTFC, LDIS]

SAC_N="Seidel Aberration Coefficents"
SCW_N="Seidel coefficients in waves"
TAC_N="Transverse Aberration Coefficents"
LAC_N="Longitudinal Aberration Coefficents"

AberNames=[SAC_N, SCW_N, TAC_N, LAC_N]
return AberNames, SAC, SCW, TAC, LAC
"""




# ############################################################

Pup = kn.pupilcalc(Doblete, sup, W, AperType, AperVal)
Pup.Samp = 25
Pup.Ptype = "fan"
Pup.FieldY = field

x, y, z, L, M, N = Pup.Pattern2Field()
Rayos = kn.raykeeper(Doblete)

for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()

kn.display2d(Doblete, Rayos, 0)