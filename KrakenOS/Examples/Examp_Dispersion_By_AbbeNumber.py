#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Doublet Lens"""


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

# ______________________________________#


# From Examp_Doublet_Lens.py
P_Obj = Kos.surf(Rc=0.0, Thickness=10, Glass='AIR', Diameter=30.0)
L1a = Kos.surf(Rc=9.284706570002484E+001, Thickness=6.0, Glass='BK7', Diameter=30.0, Axicon=0)
L1b = Kos.surf(Rc=-3.071608670000159E+001, Thickness=3.0, Glass='F2', Diameter=30.0)

L1c = Kos.surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass='AIR', Diameter=30.0)
# L1c = Kos.surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass=1, Diameter=30.0)
# L1c = Kos.surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass='nvk 1.0 0 0', Diameter=30.0)
# L1c = Kos.surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass='___BLANK 1 0 1.52216 5.88E+1 0 0 0 0 0 0'  , Diameter=30.0)

'''
I added input format for Glass attribute in surf class
Physics.py
KrakenSys.py
file are modified

case I)
as like the original version, input is material name [str]

ex)
L1c = surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass='AIR', Diameter=30.0)

case II)
as like the original version, input is just refractive index [float, int]

ex)
L1c = surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass=1, Diameter=30.0)

case III)
Input is refractive index, dispersion (abbe number) and alpha
wavelength dependency can be corrected by using Abbe number
I refered the dispersion correction equation from
https://github.com/quartiq/rayopt/blob/master/rayopt/material.py

format
'nvk <refractive index> <Abbe number> <Alpha>'

ex) 'nvk 1.0 0 0' [str]
L1c = surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass='nvk 1.0 0 0', Diameter=30.0)

case IIII)
When implement the zemax lens data, there is a BLANK material that only refractive index and dispersion are recorded
To read that format, I modified code.

reference for the zemax glass data format
https://github.com/nzhagen/zemaxglass/blob/master/ZemaxGlass_user_manual.pdf
NM <glass_name> <dispersion_formula_number> <MIL> <N(d)> <V(d)> <Exclude_sub> <status> <melf_freq>

format
'___BLANK 1 0 <refractive index> <Abbe number> 0 0 0 0 0 0'

ex) '___BLANK 1 0 1.52216 5.88E+1 0 0 0 0 0 0'
L1c = surf(Rc=-7.819730726078505E+001, Thickness=9.737604742910693E+001, Glass='___BLANK 1 0 1.52216 5.88E+1 0 0 0 0 0 0'  , Diameter=30.0)
'''
P_Ima = Kos.surf(Rc=0.0, Thickness=0.0, Glass='AIR', Diameter=50.0, Name='Plano imagen')
# ______________________________________#

# ______________________________________#

A = [P_Obj, L1a, L1b, L1c, P_Ima]
config_1 = Kos.Setup()

# ______________________________________#

Doblete = Kos.system(A, config_1)
Rayos1 = Kos.raykeeper(Doblete)
Rayos2 = Kos.raykeeper(Doblete)
Rayos3 = Kos.raykeeper(Doblete)
RayosT = Kos.raykeeper(Doblete)

# ______________________________________#

tam = 10
rad = 10.0
tsis = len(A) - 1
for j in range(-tam, tam + 1):
    for i in range(-tam, tam + 1):
        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
            W = 0.4
            Doblete.Trace(pSource_0, dCos, W)
            Rayos1.push()
            RayosT.push()
            W = 0.5
            Doblete.Trace(pSource_0, dCos, W)
            Rayos2.push()
            RayosT.push()
            W = 0.6
            Doblete.Trace(pSource_0, dCos, W)
            Rayos3.push()
            RayosT.push()

# ______________________________________#

Kos.display2d(Doblete, RayosT, 0)


