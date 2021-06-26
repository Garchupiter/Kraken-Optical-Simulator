#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 12:04:14 2020

@author: joelherreravazquez
"""
import os
import sys

import numpy as np

import Kraken as kn

currentDirectory = os.getcwd()
sys.path.insert(1, currentDirectory + '/library')

#############################################################


##############################################################    
P_Obj = kn.surf()
P_Obj.Rc = 0.0
P_Obj.Thickness = 100
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0
P_Obj.Name = "P Obj"

L1a = kn.surf()
L1a.Rc = 9.284706570002484E+001
L1a.Thickness = 6.0
L1a.Glass = "N-BK7"
L1a.Diameter = 30.0
L1a.Axicon = 0

L1b = kn.surf()
L1b.Rc = -3.071608670000159E+001

L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

L1c = kn.surf()
L1c.Rc = -7.819730726078505E+001
L1c.Thickness = 9.737604742910693E+001 - 40
L1c.Glass = "AIR"
L1c.Diameter = 30

pupila = kn.surf()
pupila.Rc = 0
pupila.Thickness = 40.
pupila.Glass = "AIR"
pupila.Diameter = 15.0
pupila.Name = "Ap Stop"

P_Ima = kn.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 20.0
P_Ima.Name = "Plano imagen"

A = [P_Obj, L1a, L1b, L1c, pupila, P_Ima]

##########################################################################


config_1 = kn.Kraken_setup()
Doblete = kn.system(A, config_1)

SYSTEM = Doblete
W = 0.6
sup = 4
AperVal = 3
AperType = "STOP"

Pup = kn.pupilcalc(SYSTEM, sup, W, AperType, AperVal)
Prx = SYSTEM.Parax(W)
SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, C, n, d = Prx

n = np.asarray(n)
n = np.insert(n, 0, 1)
n_elements = SYSTEM.n

c = []
d = []
kon = []
d.append(0)

for i in range(0, n_elements):
    rc = SYSTEM.SDT[i].Rc
    if rc == 0:
        c.append(rc)
    else:
        c.append(1 / rc)

    t = SYSTEM.SDT[i].Thickness
    d.append(t)
    kon.append(SYSTEM.SDT[i].k)

d[1] = d[1] - Pup.PosPupInp[2]

print("Radio Pupila", Pup.RadPupInp)

Diam = Pup.RadPupInp * 2.0
teta = 0.025
c = np.asarray(c)

epsilon = kon  # [-1.07731,  -4.3281, 0, 0, 0]

asp4 = [0, 0, 0, 0, 0, 0, 0]
asp2 = [0, 0, 0, 0, 0, 0]

asp2 = np.asarray(asp2)
lamb = 0.0004

##### recalculando la contribucion de la esferica en c
""" Warren Smith - Modern optical engineering_ the design of optical systems-McGraw Hill (2008) pag 112"""
c = c + 2.0 * asp2

##########################################

u = []
h = []
u_bar = []
h_bar = []
E = []

Delta_unn = []
Delta_nn = []
Delta_n = []
Delta_1sn = []

A = []
A_bar = []

sI = []
sII = []
sIII = []
sIV = []
sV = []
sVI = []

sI_k = []
sII_k = []
sIII_k = []
sIV_k = []
sV_k = []
sVI_k = []

sI_a = []
sII_a = []
sIII_a = []
sIV_a = []
sV_a = []
sVI_a = []
# ---------------------------------------

u.append(0)
h.append(Diam / 2.0)

u_bar.append(np.deg2rad(teta))
h_bar.append(0)
E.append(0)

for k in range(0, len(c) - 1):
    nk = n[k]
    nkm1 = n[k + 1]
    if np.abs(nk) != np.abs(nkm1):
        nk = np.abs(nk)

    u_pro = (nk * u[k] - (h[k] * c[k] * (nkm1 - nk))) / nkm1
    u.append(u_pro)
    h_pro = h[k] + (u[k + 1] * d[k + 1])
    h.append(h_pro)

    ub_pro = (nk * u_bar[k] - (h_bar[k] * c[k] * (nkm1 - nk))) / nkm1
    u_bar.append(ub_pro)
    hb_pro = h_bar[k] + (u_bar[k + 1] * d[k + 1])
    h_bar.append(hb_pro)

    Epro = (-d[k + 1] / (nkm1 * h[k + 1] * h[k])) + E[k]
    E.append(Epro)

for k in range(0, len(c) - 1):

    nk = n[k]
    nkm1 = n[k + 1]
    if np.abs(nk) != np.abs(nkm1):
        nk = np.abs(nk)

    UNNN = (u[k + 1] / nkm1) - (u[k] / nk)
    Delta_unn.append(UNNN)

    NN = (1.0 / nkm1) - (1.0 / nk)
    Delta_nn.append(NN)

    N = nkm1 - nk
    Delta_n.append(N)

    s1N = (1.0 / (nkm1 ** 2)) - (1 / (nk ** 2))
    Delta_1sn.append(s1N)

    A.append(nk * ((h[k] * c[k]) + u[k]))

Abar = n[0] * (h_bar[0] * c[0] + u_bar[0])
H = (A[0] * h_bar[0]) - (Abar * h[0])

for k in range(0, len(c) - 1):
    A_bar.append((H / h[k]) * ((A[k] * h[k] * E[k]) - 1))

for k in range(0, len(c) - 1):
    s1 = (A[k] ** 2.) * h[k] * Delta_unn[k]
    sI.append(s1)

    s2 = A[k] * A_bar[k] * h[k] * Delta_unn[k]
    sII.append(s2)

    s3 = (A_bar[k] ** 2.) * h[k] * Delta_unn[k]
    sIII.append(s3)

    s4 = (H ** 2.) * c[k] * Delta_nn[k]
    sIV.append(s4)

    # asdf=(A_bar[k]**2.0) * Delta_1sn[k] * h[k]
    # P = c[k] * Delta_nn[k]
    # lkj = (E[k] + (A_bar[k] * h[k]))
    # s5 = A_bar[k] * ( asdf -  (lkj * h_bar[k] * P) )  

    if A[k] != 0:

        H
        P = c[k] * Delta_nn[k]
        AbarA = A_bar[k] / A[k]

        p1 = (H ** 2.) * P
        p2 = (A_bar[k] ** 2.) * h[k] * Delta_unn[k]

        s5 = AbarA * (p1 + p2)
    else:
        s5 = 0

    sV.append(s5)

    # s5_0 = ((()  * h[k]) * Delta_unn[k]
    # s5_1 = ((A_bar[k] / A[k])  * H * H * c[k]) * Delta_nn[k]
    # s5= s5_0
    # sV.append(s5)

    # ------------------------------------------------
    a = epsilon[k] * (c[k] ** 3.) * (h[k] ** 4.)

    s1k = a * Delta_n[k] * ((h_bar[k] / h[k]) ** 0)
    sI_k.append(s1k)

    s2k = a * Delta_n[k] * ((h_bar[k] / h[k]) ** 1)
    sII_k.append(s2k)

    s3k = a * Delta_n[k] * ((h_bar[k] / h[k]) ** 2)
    sIII_k.append(s3k)

    s3k = a * 0.
    sIII_k.append(s3k)

    s5k = a * Delta_n[k] * ((h_bar[k] / h[k]) ** 3)
    sIV_k.append(s5k)

    s6k = a * Delta_n[k] * ((h_bar[k] / h[k]) ** 4)
    sIV_k.append(s6k)

    # ------------------------------------------------

    # s1a = aspher[k] * (8.) * (h[k]**4.) * Delta_n[k]
    # sI_a.append(s1a)

    # s2a = aspher[k] * (8.) * (h[k]**4.) * Delta_n[k] * (h_bar[k] / h[k])
    # sII_a.append(s2a)

    # s3a = aspher[k] * (8.) * (h[k]**4.) * Delta_n[k] * ((h_bar[k] / h[k])**2.)
    # sIII_a.append(s3a)

    Asp = asp4[k] - ((asp2[k] / 4.0) * ((4.0 * asp2[k] ** 2) + (6.0 * asp2[k] * c[k]) + (3.0 * c[k] ** 2.0)))

    s1a = Asp * 8. * (h[k] ** 4.) * Delta_n[k]
    sI_a.append(s1a)

    s2a = Asp * 8. * (h[k] ** 4.) * Delta_n[k] * (h_bar[k] / h[k])
    sII_a.append(s2a)

    s3a = Asp * 8. * (h[k] ** 4.) * Delta_n[k] * (h_bar[k] / h[k] ** 2.)
    sIII_a.append(s3a)

sI = np.asarray(sI)
sII = np.asarray(sII)
sIII = np.asarray(sIII)
sIV = np.asarray(sIV)
sV = np.asarray(sV)

sI_k = np.asarray(sI_k)
sII_k = np.asarray(sII_k)
sIII_k = np.asarray(sIII_k)
sIV_k = np.asarray(sIV_k)

SI = (np.sum(sI) - np.sum(sI_k) - np.sum(sI_a)) / (8. * lamb)
SII = (np.sum(sII) - np.sum(sII_k) - np.sum(sII_a)) / (2. * lamb)
SIII = (np.sum(sIII) - np.sum(sIII_k) - np.sum(sIII_a)) / (2. * lamb)
SIV = (np.sum(sIV)) / (4. * lamb)
SV = (np.sum(sV) - np.sum(sIV_k)) / (2. * lamb)

print(SI, SII, SIII, SIV, SV)

AperVal = 3
AperType = "EPD"
Pup = kn.pupilcalc(Doblete, sup, W, AperType, AperVal)

Pup.Samp = 25
Pup.Ptype = "fan"
Pup.FieldType = "angle"
Pup.FieldY = 1.0

x, y, z, L, M, N = Pup.Pattern2Field()
Rayos = kn.raykeeper(Doblete)

for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]

    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()

# X,Y,Z,L,M,N=Rayos.pick(-1)
# plt.figure(300)
# plt.plot(X,Y, 'x')
# plt.axis('square')
# plt.show(block = False)


Pup.FieldY = -1.0
x, y, z, L, M, N = Pup.Pattern2Field()
for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]

    Doblete.Trace(pSource_0, dCos, W)
    Rayos.push()

kn.display2d(Doblete, Rayos, 0)

# plt.figure(300)
# plt.plot(x,y, 'x')
# plt.axis('square')
# plt.show(block = False)
