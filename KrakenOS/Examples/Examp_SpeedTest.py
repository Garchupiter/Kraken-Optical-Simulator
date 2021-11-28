#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens Pupil Seidel"""
from numba import jit
import timeit
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

# _________________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 100
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0
P_Obj.Name = "P Obj"

# _________________________________________#

L1a = Kos.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "N-BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# _________________________________________#

L1b = Kos.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# _________________________________________#

L1c = Kos.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001 - 40
L1c.Glass = "AIR"
L1c.Diameter = 30

# _________________________________________#

pupila = Kos.surf()
pupila.Rc = 0
pupila.Thickness = 40.
pupila.Glass = "AIR"
pupila.Diameter = 15.0
pupila.Name = "Ap Stop"

# _________________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 20.0

# _________________________________________#

A = [P_Obj, L1a, L1b, L1c, pupila, P_Ima]
config_1 = Kos.Setup()

# _________________________________________#




Doblete = Kos.system(A, config_1)

# _________________________________________#

W = 0.6
Surf = 4
AperVal = 3
AperType = "EPD"
fieldType = "angle"

Pup = Kos.PupilCalc(Doblete, Surf, W, AperType, AperVal)

Pup.Samp = 45
Pup.Ptype = "hexapolar"
Pup.FieldY = 3.25


#_________________________________________#

x, y, z, L, M, N = Pup.Pattern2Field()
Rayos = Kos.raykeeper(Doblete)

# _________________________________________#



start = timeit.default_timer()
print("-----------------")
cont=0.0
for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()
    cont = cont + 1.0

    # print(Doblete.NP_XYZ[0:Doblete.iter+1])
    # print(Doblete.NP_S_XYZ[0:Doblete.iter+1])




stop = timeit.default_timer()
execution_time = stop - start

print("Program Executed in "+str(execution_time/(cont*5.0))+" per ray segment")
print(cont*5, " Rays segments trace d")


# 25-oct-21 Program Executed in 0.0004910416639832562 31055.0  Segmentos
# 27-oct-21 Program Executed in 0.0002906117962646915 31055.0  Segmentos
# 29-oct-21 Program Executed in 0.00018554626739655485 31055.0  Segmentos
# 29-oct-21 Program Executed in 0.0001726453898889197 31055.0  Segmentos
# 29-oct-21 Program Executed in 0.00016649437021413582 31055.0  Segmentos
# 29-oct-21 Program Executed in 0.00015742018998550996 31055.0  Segmentos
# 30-oct-21 Program Executed in 0.0001484923241024164 31055.0  Segmentos

# 31-oct-21 implementación de FastTrace
# 31-oct-21 Program Executed in 0.00013013021978017754 455.0  Segmentos
# 31-oct-21 Program Executed in 0.00011787373611336453 31055.0  Segmentos
# 31-oct-21 Program Executed in 0.00011517481751730815 31055.0  Segmentos

# cambio de np.linalg.norm a np.sqrt(sum x**2)
# 01-nov-21 Program Executed in 0.000108276663693441 31055.0  Segmentos


# _________________________________________#
# Kos.display3d(Doblete, Rayos, 0)
