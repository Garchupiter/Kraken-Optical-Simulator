
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

from scipy import optimize
from scipy.optimize import fmin_cg



import Kraken as kn

import numpy as np
import matplotlib.pyplot as plt
import time




##############################################################    
##############################################################    
P_Obj=kn.surf()
P_Obj.Rc=0
P_Obj.Thickness=1000+3.452200000000000E+003
P_Obj.Glass="AIR"
P_Obj.Diameter=1.059E+003*2.0

Thickness=3.452200000000000E+003
M1=kn.surf()
M1.Rc=-9.638000000004009E+003
M1.Thickness=-Thickness
M1.k=-1.077310000000000E+000
M1.Glass="MIRROR"
M1.Diameter=1.059E+003*2.0
M1.InDiameter=250*2.0

M2=kn.surf()
M2.Rc=-3.93E+003
M2.Thickness=Thickness+1.037525880125084E+003
M2.k=-4.328100000000000E+000
M2.Glass="MIRROR"
M2.Diameter=3.365E+002*2.0
M2.AxisMove=0

P_Ima=kn.surf()
P_Ima.Diameter=600.0
P_Ima.Glass="AIR"
P_Ima.Name="Plano imagen"

A=[P_Obj,M1,M2,P_Ima]
configuracion_1=kn.Kraken_setup()



Telescopio=kn.system(A,configuracion_1)

Rayos0 = kn.raykeeper(Telescopio)



sup=1
W=0.55
Telescopio.Parax(W)

print(Telescopio.EFFL)
Pup=kn.pupilcalc(Telescopio, sup, W)
Telescopio.Vignetting(0)

#######################################

Pup.Samp=1
Pup.Ptype = "fany" 
Pup.FieldType = "angle"
Pup.FieldX=0.
rays=Pup.Pattern2Field()


# Pup.FieldY=0.141
# raysP4=Pup.Pattern2Field()
# raysP4=np.copy(raysP4)

# Pup.FieldY=0.0
# raysC0=Pup.Pattern2Field()
# raysC0=np.copy(raysC0)

# Pup.FieldY=0.0997
# raysM4=Pup.Pattern2Field()
# raysM4=np.copy(raysM4)





# ############################################




def TAOS_ZCP(AAA, BBB, r):
    W=0.4

    x,y,z,l,m,n=r
    respM2_DespY = np.copy(M2.DespY)
    respM2_TiltX = np.copy(M2.TiltX)
    
    M2.DespY = BBB
    M2.TiltX = AAA
    Telescopio.ResetData()
    Telescopio.IgnoreVignetting(0)
    
    for i in range(0,len(x)):
        pSource_0 = [x[i], y[i], z[i]]
        dCos=[l[i], m[i], n[i]]
        Telescopio.Trace(pSource_0, dCos, W)
        
        Rayos0.push()
    
    X,Y,Z,L,M,N=Rayos0.pick(-1)
    
       
    Ymb = Y[0]
    Zmb = Z[0]
    Mmb = M[0]
    Nmb = N[0]
    
    Yma = Y[2]
    Zma = Z[2]
    Mma = M[2]
    Nma = N[2]
    
    Zp0 = Z[1]
    Yp0 = Y[1]
    Mp = M[1]
    Np = N[1]
    
    Ym=((Mma*Mmb)/(Nma*Mmb-Nmb*Mma)    )*(Zmb-Zma+((Yma*Nma/Mma)-(Ymb*Nmb/Mmb)))
    Zm=((Ym-Ymb)*Nmb/Mmb)+Zmb
    
    Yp = Yp0 + (Mp / Np * (Zm - Zp0))
    delta_y = np.abs(Yp - Ym)
    
    # print("asdfgh : ",respM2_DespY, M2.DespY , respM2_TiltX, M2.TiltX )
    M2.DespY = respM2_DespY
    M2.TiltX = respM2_TiltX

    
    Telescopio.ResetData()
    Telescopio.IgnoreVignetting(0)
    Rayos0.clean()

    return delta_y





def Newthon_Raphson(M2DY,r):
    h=0.0001
    
    
    M2TX_1=0.0
    while True:
        Der_TAOS=(TAOS_ZCP(M2TX_1+h, M2DY, r)-TAOS_ZCP(M2TX_1-h, M2DY, r))/(2*h)
        F=TAOS_ZCP(M2TX_1, M2DY, r)
        M2TX_2=M2TX_1-(F/Der_TAOS)
        
        
        if np.abs(F)<0.0001:
            break
        M2TX_1=M2TX_2
        
       
    return M2TX_1



for i in range(1,5):
    
    M2DY=i    
    M2TX=Newthon_Raphson(M2DY,rays)
    
    ZCP=M2DY/np.sin(np.deg2rad(M2TX))
    print(M2DY, M2TX, ZCP)


M2.DespY = M2DY
M2.TiltX = M2TX

# M2.TiltX=0.4519959691826365
# M2.DespY=3.0

Telescopio.ResetData()
Telescopio.ResetSolid()  
Telescopio.IgnoreVignetting(0)


# ###########################################



x,y,z,l,m,n=rays

for i in range(0,len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos=[l[i], m[i], n[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    
    Rayos0.push()





# x,y,z,l,m,n=raysC0

# for i in range(0,len(x)):
#     pSource_0 = [x[i], y[i], z[i]]
#     dCos=[l[i], m[i], n[i]]
#     Telescopio.Trace(pSource_0, dCos, W)
    
#     Rayos0.push()





# x,y,z,l,m,n=raysM4

# for i in range(0,len(x)):
#     pSource_0 = [x[i], y[i], z[i]]
#     dCos=[l[i], m[i], n[i]]
#     Telescopio.Trace(pSource_0, dCos, W)
    
#     Rayos0.push()


# ############################################
      
kn.display2d(Telescopio,Rayos0,0)


# X,Y,Z,L,M,N=Rayos0.pick(-1)

# cenX=np.mean(X)
# cenY=np.mean(Y)

# plt.plot(X,Y, 'x')



# # axis labeling
# plt.xlabel('numbers')
# plt.ylabel('values')

# # figure name
# plt.title('Dot Plot : Red Dots')
# plt.axis('square')
# plt.show()   


