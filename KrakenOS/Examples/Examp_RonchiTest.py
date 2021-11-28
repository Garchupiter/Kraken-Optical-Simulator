import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import scipy


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
######################################################################################


P_Obj = Kos.surf()

P_Obj.Thickness = 1.0
P_Obj.Glass = "AIR"
P_Obj.Diameter = 100
P_Obj.Drawing = 0



Ronchi = Kos.surf()
Ronchi.Diameter = 100
Ronchi.Thickness = 5000
Ronchi.Glass = "AIR"
Sp=0.5
n= np.arange(-10,10,Sp)
AAA = pv.MultiBlock()

for j in n:
    Rulling = pv.Plane(center=[j, 0, 0], direction=[0, 0, 1], i_size=Sp/2, j_size=110, i_resolution=5, j_resolution=5)
    AAA.append(Rulling)
Ronchi.Mask_Shape = AAA
Ronchi.Mask_Type = 2


Mirror = Kos.surf()
Mirror.Rc = -5000
Mirror.k = -0.03
Mirror.Diameter = 2500
Mirror.Thickness = -5000
Mirror.Glass = "MIRROR"



Ronchi2 = Kos.surf()
Ronchi2.Diameter = 100
Ronchi2.Thickness = -5000
Ronchi2.Glass = "AIR"
Ronchi2.Mask_Shape = AAA
Ronchi2.Mask_Type = 2



P_Ima = Kos.surf()
P_Ima.Rc = 0
P_Ima.Thickness = -.0001
P_Ima.Glass = "AIR"
P_Ima.Diameter = 6000.0
P_Ima.Drawing = 1
P_Ima.Name = "Plano imagen"

A = [P_Obj, Ronchi, Mirror, Ronchi2, P_Ima]

configur = Kos.Setup()
Telescope = Kos.system(A, configur)
Rays = Kos.raykeeper(Telescope)

W = 0.633
Sun = Kos.SourceRnd()

def f(x):
    res = 1
    return res


Sun.field = 15.0

Sun.fun = f
Sun.dim = 1.5
Sun.num = 10000
L, M, N, X, Y, Z = Sun.rays()

Xr = np.zeros_like(L)
Yr = np.zeros_like(L)
Nr = np.zeros_like(L)

con = 0
con2 = 0
for i in range(0, Sun.num):
    if con2 == 10:
        print(100. * i / Sun.num)
        con2 = 0

    pSource_0 = [X[i], Y[i], Z[i]]
    dCos = [L[i], M[i], N[i]]
    Telescope.Trace(pSource_0, dCos, W)

    Xr[con] = Telescope.Hit_x[-1]
    Yr[con] = Telescope.Hit_y[-1]
    Nr[con] = Telescope.SURFACE[-1]
    con = con + 1
    con2 = con2 + 1
    # Rays.push()

args = np.argwhere(Nr == len(A)-1)
plt.plot(Xr[args], Yr[args], '.', c="g", markersize=1)

# axis labeling
plt.xlabel('x')
plt.ylabel('y')

# figure name
plt.title('Spot diagram')
plt.axis('square')
plt.show()

#             Rays.push()
# Kos.display3d(Telescope, Rays, 1)
