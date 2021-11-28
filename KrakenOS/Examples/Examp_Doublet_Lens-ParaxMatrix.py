#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens Para xMatrix"""

import time
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

start_time = time.time()

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 10
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# ______________________________________#

L1a = Kos.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

# ______________________________________#

L1b = Kos.surf()
L1b.Rc = -3.071608670000159E+001
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# ______________________________________#

L1c = Kos.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001
L1c.Glass = "AIR"
L1c.Diameter = 30

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 3.0
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, L1a, L1b, L1c, P_Ima]
config_1 = Kos.Setup()

# ______________________________________#

Doblete = Kos.system(A, config_1)
Prx = Doblete.Parax(0.4)
SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, CC, N_Prec, DD = Prx
print(EFFL)

# ______________________________________#

L1a.Rc = L1a.Rc + 1
Doblete.SetData()
Prx = Doblete.Parax(0.4)
SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, CC, N_Prec, DD = Prx
print(EFFL)

# ______________________________________#

L1a.Rc = L1a.Rc + 1
Doblete.SetData()
Prx = Doblete.Parax(0.4)
SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, CC, N_Prec, DD = Prx
print(EFFL)

# ______________________________________#


L1a.Rc = L1a.Rc + 1
Doblete.SetData()
Prx = Doblete.Parax(0.4)
SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, CC, N_Prec, DD = Prx
print(EFFL)

print("=======================================")

print(SistemMatrix)
print(S_Matrix)
print( N_Matrix)
print( a)
print( b )
print( c)
print( d)
print( EFFL)
print( PPA)
print( PPP)
print( CC)
print( N_Prec)
print( DD)
