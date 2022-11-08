#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Diffraction Grating Transmission"""

import numpy as np
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
# _________________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 5
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# _________________________________________#

Prism_c1 = Kos.surf()
Prism_c1.Rc = 0.0
Prism_c1.Thickness = 13
Prism_c1.Glass = "BK7"
Prism_c1.Diameter = 30.0
Prism_c1.TiltX=20
Prism_c1.AxisMove = 0

# _________________________________________#

Prism_c2 = Kos.surf()
Prism_c2.Rc = 0.0
Prism_c2.Thickness = 30
Prism_c2.Glass = "AIR"
Prism_c2.Diameter = 30.0
Prism_c2.TiltX = -20
Prism_c2.AxisMove = 0
# _________________________________________#


P_Ima = Kos.surf()
P_Ima.Name = "Plano imagen"
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 100.0

# _________________________________________#

A = [P_Obj, Prism_c1, Prism_c2,  P_Ima]
configuracion_1 = Kos.Setup()

# _________________________________________#

Doblete = Kos.system(A, configuracion_1)


for ii in range(0,40):
    Rayos = Kos.raykeeper(Doblete)

    ii = 12
    tet = ii
    x_0 = 0
    y_0 = 0

    pSource_0 = [x_0, y_0, 0.0]
    dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
    W = 0.4
    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()
    W = 0.5
    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()
    W = 0.6
    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()

    # _________________________________________#



    x0,y0,z0,l0,m0,n0 = Rayos.pick(1)
    x1,y1,z1,l1,m1,n1 = Rayos.pick(-1,coordinates="local")

    S0 = np.sqrt((l0**2) + (m0**2) + (n0**2))
    S1 = np.sqrt((l1**2) + (m1**2) + (n1**2))
    CosSigma = ((l0*l1) + (m0*m1) + (n0*n1))/ (S0*S1)
    sigma = np.arccos(CosSigma)


    alpha = np.deg2rad(40)


    print("Sigma: ", ii, np.rad2deg(sigma))




N = np.sin((1/2) * (sigma + alpha)) / np.sin((1/2) * alpha)
print(N)



[SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, APP, PPP, C, N, D] = Doblete.Parax(0.4)
print(N)
[SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, APP, PPP, C, N, D] = Doblete.Parax(0.5)
print(N)
[SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, APP, PPP, C, N, D] = Doblete.Parax(0.6)
print(N)



Kos.display2d(Doblete, Rayos, 0)
