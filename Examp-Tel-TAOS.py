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


P_Obj.Thickness= 2.000000000000000E+003
P_Obj.Glass= "AIR"
P_Obj.Diameter=(6.796727741707513E+002*2.0)
P_Obj.Drawing= 0


M1=kn.surf()
M1.Rc=(-6.06044E+003)
M1.Thickness=(-1.774190000000000E+003+1.853722901194000E+000)
M1.k=(-1.637E+000)
M1.Glass= "MIRROR"
M1.Diameter=(6.63448E+002*2.0)
M1.InDiameter=(228.6*2.0)
M1.DespY= 0.0
M1.TiltX= 0.0
M1.AxisMove= 0

M2=kn.surf()
M2.Rc=(-6.06044E+003)
M2.Thickness=(-M1.Thickness)
M2.k=(-3.5782E+001)
M2.Glass= "MIRROR"
M2.Diameter=(2.995730651164167E+002*2.0)

ED0=np.zeros(20)

ED0[2]=4.458178314555000E-018
M2.ExtraData= ED0


Vertex=kn.surf()
Vertex.Thickness= 30.0
Vertex.Glass= "AIR"
Vertex.Diameter= 600.0
Vertex.Drawing= 0

Corrector_c1=kn.surf()
Corrector_c1.Thickness= 6.6E+000
Corrector_c1.Glass= "SILICASCHOTT"  #"BK7"#"LITHOSIL-Q"
Corrector_c1.Diameter=(118.0*2)
ED1=np.zeros(20)
ED1[0]=2.059727552003000E-005
ED1[1]=-1.135080732384000E-009
Corrector_c1.ExtraData= ED1


Corrector_c2=kn.surf()
Corrector_c2.Thickness=(341.65484183207997-100.0)
Corrector_c2.Glass= "AIR"
Corrector_c2.Diameter=(118.0*2.0)


Corrector_c2=kn.surf()
Corrector_c2.Thickness= 341.65484183207997
Corrector_c2.Glass= "AIR"
Corrector_c2.Diameter=(118.0*2)

P_Ima=kn.surf()
P_Ima.Rc= 0.0
P_Ima.Thickness= 0.0
P_Ima.Glass == "AIR"
P_Ima.Diameter= 300.0
P_Ima.Drawing= 0
P_Ima.Name= "Plano imagen"

############################################
configuracion_1=kn.Kraken_setup()

A=[P_Obj,M1,M2,Vertex,Corrector_c1, Corrector_c2 ,P_Ima]
Telescopio=kn.system(A,configuracion_1)
Rayos=kn.raykeeper(Telescopio)

############################################

sup=1
W=0.4
Pup=kn.pupilcalc(Telescopio, sup, W)

Pup.Samp=10
Pup.Ptype = "square" 
Pup.FieldType = "angle"

Pup.FieldX=0.0
xa,ya,za,La,Ma,Na=Pup.Pattern2Field()


############################################    


for i in range(0,len(xa)):
    pSource_0 = [xa[i], ya[i], za[i]]
    dCos=[La[i], Ma[i], Na[i]]
    Telescopio.Trace(pSource_0, dCos, W)
    Rayos.push()
    
  
    
############################################
        
kn.display3d(Telescopio,Rayos,1)

############################################



X,Y,Z,L,M,N=Rayos.pick(-1)

cenX=np.mean(X)
cenY=np.mean(Y)

plt.plot(X,Y, '.')



# axis labeling
plt.xlabel('numbers')
plt.ylabel('values')

# figure name
plt.title('Dot Plot : Red Dots')
plt.axis('square')
plt.show()   








