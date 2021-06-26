
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""

import Kraken as kn
import numpy as np
import matplotlib.pyplot as plt
import time
 
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

######################



configuracion_1=kn.Kraken_setup()

Telescopio=kn.system(A,configuracion_1)
Rayos=kn.raykeeper(Telescopio)

############################################

sup=1
W=0.4
Pup=kn.pupilcalc(Telescopio, sup, W)

Pup.Samp=7
Pup.Ptype = "fan" 
Pup.FieldType = "angle"

Pup.FieldX=0.0
xa,ya,za,La,Ma,Na=Pup.Pattern2Field()

Pup.FieldX=1.0
xb,yb,zb,Lb,Mb,Nb=Pup.Pattern2Field()

Pup.FieldX=-1.0
xc,yc,zc,Lc,Mc,Nc=Pup.Pattern2Field()

############################################    


for i in range(0,len(xa)):
    pSource_0 = [xa[i], ya[i], za[i]]
    dCos=[La[i], Ma[i], Na[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()
    
    pSource_0 = [xb[i], yb[i], zb[i]]
    dCos=[Lb[i], Mb[i], Nb[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()
    
    pSource_0 = [xc[i], yc[i], zc[i]]
    dCos=[Lc[i], Mc[i], Nc[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()
  
    
############################################        
        
kn.display3d(Telescopio,Rayos,2)



X,Y,Z,L,M,N=Rayos.pick(-1)

plt.plot(X,Y, 'x')



# axis labeling
plt.xlabel('numbers')
plt.ylabel('values')

# figure name
plt.title('Dot Plot : Red Dots')
plt.axis('square')
plt.show()        