# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Examp Tel 2M Wavefront Fitting"""

import os
import sys
import matplotlib.pyplot as plt
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

currentDirectory = os.getcwd()
sys.path.insert(1, currentDirectory + '/library')

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 1000*0 + 3452.2
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1059. * 2.0
P_Obj.Drawing=0

# ______________________________________#

Thickness = 3452.2
M1 = Kos.surf()
M1.Rc = -9638.0
M1.Thickness = -Thickness
M1.k = -1.07731
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
#M1.InDiameter = 250 * 2.0

# ______________________________________#

M2 = Kos.surf()
M2.Rc = -3930.0
M2.Thickness = 4489.731761 -1.107e-9 # Thickness + 1037.525880
M2.k = -4.3281
M2.Glass = "MIRROR"
M2.Diameter = 336.5 * 2.0


# ______________________________________#


N1 = Kos.surf()
N1.Thickness = 720
N1.Glass = "NULL"
N1.Diameter = 100
N1.TiltX = 10
N1.Order = 0
N1.AxisMove = 1

############################

Colim = Kos.surf()
Colim.Rc = -1440.0
Colim.k = -1.0
Colim.Thickness = -720
Colim.Glass = "MIRROR"
Colim.Diameter = 400


############################
N2 = Kos.surf()
N2.Glass = "NULL"
N2.DespY =125.9836755
N2.TiltX =3.181E-15
N2.AxisMove = 1
############################

N3 = Kos.surf()
N3.Glass = "NULL"
N3.TiltX =-71
N3.AxisMove = 1

#############################

G1 = Kos.surf()
G1.Glass = "MIRROR"
G1.Diameter = 300
G1.Grating_D = (1/0.079)
G1.Diff_Ord = -41
G1.Grating_Angle = 0.
G1.AxisMove=1
#############################

N4 = Kos.surf()
N4.Glass = "NULL"
N4.TiltX = -N3.TiltX
N4.AxisMove = 1

#############################

N5 = Kos.surf()
N5.Glass = "NULL"
N5.Thickness = 550
N5.TiltX = -N1.TiltX
N5.AxisMove = 1

#############################

N6 = Kos.surf()
N6.Glass = "NULL"
N6.TiltX = -4.95736765
N6.DespY = 47.70641809
N6.AxisMove = 1

#############################

N7 = Kos.surf()
N7.Glass = "NULL"
N7.TiltZ = 90
N7.AxisMove = 1
N7.Drawing = 1
N7.Diameter = 100

#############################

N8 = Kos.surf()
N8.Glass = "NULL"
N8.TiltX = -20.7
N8.AxisMove = 1

#############################

G2 = Kos.surf()
G2.Thickness = -100
G2.Glass = "MIRROR"
G2.Diameter = 300
G2.Grating_D =1/0.3
G2.Diff_Ord = 1
G2.Grating_Angle = 0
G2.SubAperture = [1.0,0,0]

#############################

N9 = Kos.surf()
N9.Glass = "NULL"
N9.DespY = -60.53784
N9.TiltX = -31.1898
N9.AxisMove = 1

#############################

N10 = Kos.surf()
N10.Glass = "NULL"
N10.TiltZ = 180
N10.AxisMove = 1

#############################

Asp1 = Kos.surf()
Asp1.Glass = "F_SILICA"
Asp1.Diameter = 100
Asp1.Rc = -14925.373
Asp1.Thickness = -10

ED1=np.zeros(10)
ED1[0] = 0.0    # 2nd Order
ED1[1] = 5.1E-9  # 4nd Order
ED1[2] = 5.1E-14 # 6nd Order
Asp1.AspherData = ED1

#############################

Asp2 = Kos.surf()
Asp2.Glass = "AIR"
Asp2.Diameter = 100
Asp2.Thickness = -175

#############################

N11 = Kos.surf()
N11.Glass = "NULL"
N11.TiltX = -22.5
N11.AxisMove = 1

#############################


Flat = Kos.surf()
Flat.Thickness = 0
Flat.Glass = "MIRROR"
Flat.Diameter = 92.5 * 2

#############################

N12 = Kos.surf()
N12.Glass = "NULL"
N12.Thickness = 175
N12.TiltX = -22.5
N12.AxisMove = 1

#############################

M3 = Kos.surf()
M3.Thickness = 0
M3.Glass = "MIRROR"
M3.Diameter = 95.23 * 2
M3.Rc = -476.0

#############################

N13 = Kos.surf()
N13.Glass = "NULL"
N13.Thickness = -175 - 38
N13.AxisMove = 1

##############################

L1a = Kos.surf()
L1a.Glass = "F_SILICA"
L1a.Diameter = 25.4*2
L1a.Thickness = -7
L1a.Rc = -98.5

L1b = Kos.surf()
L1b.Glass = "AIR"
L1b.Diameter = 25.4*2
L1b.Thickness = -8

##############################

L2a = Kos.surf()
L2a.Glass = "F_SILICA"
L2a.Diameter = 25.4*2
L2a.Thickness = -4

L2b = Kos.surf()
L2b.Glass = "AIR"
L2b.Diameter = 25.4*2
L2b.Thickness = -5.40475
#############################


P_Ima = Kos.surf()
P_Ima.Diameter = 50
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, M1, M2, N1, Colim, N2, N3, G1,N4, N5, N6, N7, N8, G2, N9, N10, Asp1, Asp2, N11, Flat, N12, M3, N13, L1a, L1b, L2a, L2b, P_Ima]
configuracion_1 = Kos.Setup()
Telescope = Kos.system(A, configuracion_1)

Rays = Kos.raykeeper(Telescope)

# ______________________________________#

tam = 5
rad = 2100 / 2
tsis = len(A) - 1

# ______________________________________#

a=np.loadtxt("thar_uves.dat.txt")
n=a[:,0]
lam=a[:,1]/10000.0



a=np.arange(-35, -57, -1)

x=[]
y=[]
z=[]
for aa in a:
    Telescope.SDT[7].Diff_Ord=aa
    Telescope.SetData()
    print("Order: ",aa)

    # ______________________________________#
    for q in range(0,len(n)):
        W = lam[q]
        i=0
        j=0
        x_0 = (i / tam) * rad
        y_0 = (j / tam) * rad
        r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
        if r < rad:
            tet = 0.0
            pSource_0 = [x_0, y_0, 0.0]
            dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]

            Telescope.Trace(pSource_0, dCos, W)
            if Telescope.SURFACE[-1] == Telescope.n-1:
                xyz = Telescope.OST_XYZ[-1]
                x.append(xyz[0])
                y.append(xyz[1])
                z.append(xyz[2])
                Rays.push()

X=np.asarray(x)
Y=np.asarray(y)
Z=np.asarray(z)

# ______________________________________#


plt.plot(X, -Y, 'x')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Spot Diagram')
plt.axis('square')
plt.show()



Kos.display3d(Telescope, Rays, 0)
