import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

import Kraken as kn

P_Obj = kn.surf()

P_Obj.Thickness = 5000.0
P_Obj.Glass = "AIR"
P_Obj.Diameter = 6.796727741707513E+002 * 2.0
P_Obj.Drawing = 0

######################################################################################


objeto = kn.surf()
objeto.Rc=-12000
objeto.k=-1
objeto.Diameter=2500
objeto.Thickness = -6000
objeto.Glass = "MIRROR"


P_Ima = kn.surf()
P_Ima.Rc = 0
P_Ima.Thickness = -1.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 6000.0
P_Ima.Drawing = 1
P_Ima.Name = "Plano imagen"

A = [P_Obj, objeto, P_Ima]

configur = kn.Kraken_setup()
Telescope = kn.system(A, configur)
Rays = kn.raykeeper(Telescope)

W = 0.633


def f(x):
    Wh=0.025
    r=(x*90.0/0.025)*np.pi
    res=np.sin(r)/r
    return res

# def f(x):
#     res=1
#     return res

# def f(x):
#     r=(x*90.0/0.025)
#     res=r**2
#     return res


Sun=kn.SourceRnd()
Sun.fun=f
Sun.dim=3000
Sun.field=0.025*3
Sun.num=100000
L, M, N, X, Y, Z = Sun.rays()

         
Xr=np.zeros_like(L)
Yr=np.zeros_like(L)
Nr=np.zeros_like(L)

con=0
con2=0
for i in range(0,Sun.num):
    if con2==10:
        print(100.*i/Sun.num)
        con2=0
    
    pSource_0 = [X[i], Y[i], Z[i]]
    dCos = [L[i], M[i], N[i]]
    Telescope.Trace(pSource_0, dCos, W)

    Xr[con]=Telescope.Hit_x[-1]
    Yr[con]=Telescope.Hit_y[-1]
    Nr[con]=Telescope.SURFACE[-1]
    con=con+1
    con2=con2+1
    #Rays.push()

args=np.argwhere(Nr==2)
plt.plot(Xr[args], Yr[args], '.', c="g", markersize=1)

# axis labeling
plt.xlabel('x')
plt.ylabel('y')

# figure name
plt.title('Dot Plot')
plt.axis('square')
plt.show()




#             Rays.push()
#kn.display3d(Telescope, Rays, 0)
