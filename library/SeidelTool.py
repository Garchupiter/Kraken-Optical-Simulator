#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 08:04:35 2021

@author: joelherreravazquez
Welford, W.T. (1986). Aberrations of Optical Systems (1st ed.). Routledge. https://doi.org/10.1201/9781315136530


Warren Smith - Modern optical engineering_ the design of optical systems-McGraw Hill (2008).pdf


Zemax - Optical design program - User Manual-Juy 8 2011
https://neurophysics.ucsd.edu/Manuals/Zemax/ZemaxManual.pdf
"""

import numpy as np
import Kraken as kn
from pupil_tool import pupilcalc, SolveVectCross


def Seidel(SYSTEM, sup, W, ApType, ApVal, field, fieldType):
    Pup = pupilcalc(SYSTEM, sup, W, ApType, ApVal)
    Prx = SYSTEM.Parax(W)
    SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, C, NN, d = Prx

    n_1 = np.asarray(NN)  # prev
    n_sign = np.sign(NN)
    n_abs = np.abs(NN)

    for si in range(0, len(n_sign)):
        n_abs[si::] = n_sign[si] * n_abs[si::]

    n_1 = np.copy(n_abs)
    n_2 = np.copy(n_abs)
    n_1 = np.insert(n_1, 0, n_1[0])

    print(n_2)
    print(n_1)
    n_elements = SYSTEM.n
    c = []
    d = []
    kon = []
    asp2 = []
    asp4 = []
    asp6 = []

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

        asp2.append(SYSTEM.SDT[i].AspherData[0])
        asp4.append(SYSTEM.SDT[i].AspherData[1])
        asp6.append(SYSTEM.SDT[i].AspherData[2])

    c.append(0)
    d.append(0)
    kon.append(0)
    asp2.append(0)
    asp4.append(0)
    asp6.append(0)

    t = d[1]
    d[1] = d[1] - Pup.PosPupInp[2]
    Diam = Pup.RadPupInp * 2.0

    c = np.asarray(c)
    epsilon = np.asarray(kon)
    asp2 = np.asarray(asp2)
    asp4 = np.asarray(asp4)
    asp6 = np.asarray(asp6)

    """ ---------------- Standard and conic --------------------------------    
        Seidel sums for standard surfaces with conicity and Asp2o   
    ------------------------------------------------------------------------"""
    c_orig = np.copy(c)
    c = c + 2.0 * asp2

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

    sI_k = []
    sII_k = []
    sIII_k = []
    sIV_k = []
    sV_k = []

    # fieldType= angle, height

    if fieldType == "angle":
        teta_marg = 0
        teta_bar = -np.deg2rad(field)

    PosObj = t + d[1]
    if fieldType == "height":
        teta_marg = np.arctan2(Pup.RadPupInp, PosObj)
        teta_bar = -np.arctan2(field, PosObj)

    u.append(teta_marg)
    h.append(Diam / 2.0)

    u_bar.append(teta_bar)
    h_bar.append(0)
    E.append(0)

    for k in range(0, len(c) - 1):
        nk = n_1[k]
        nkm1 = n_2[k]

        u.append((nk * u[k] - (h[k] * c[k] * (nkm1 - nk))) / nkm1)
        h.append(h[k] + (u[k + 1] * d[k + 1]))

        u_bar.append((nk * u_bar[k] - (h_bar[k] * c[k] * (nkm1 - nk))) / nkm1)
        h_bar.append(h_bar[k] + (u_bar[k + 1] * d[k + 1]))

        Epro = (-d[k + 1] / (nkm1 * h[k + 1] * h[k])) + E[k]
        E.append(Epro)

        Delta_unn.append((u[k + 1] / nkm1) - (u[k] / nk))
        Delta_nn.append((1.0 / nkm1) - (1.0 / nk))
        Delta_n.append(nkm1 - nk)
        Delta_1sn.append((1.0 / (nkm1 ** 2)) - (1 / (nk ** 2)))

        A.append(nk * ((h[k] * c[k]) + u[k]))

    Abar = 1.0 * (h_bar[0] * c[0] + u_bar[0])

    H = (A[0] * h_bar[0]) - (Abar * h[0])

    for k in range(0, len(c) - 1):
        A_bar.append((H / h[k]) * ((A[k] * h[k] * E[k]) - 1))

        sI.append((A[k] ** 2.) * h[k] * Delta_unn[k])

        sII.append(A[k] * A_bar[k] * h[k] * Delta_unn[k])

        sIII.append((A_bar[k] ** 2.) * h[k] * Delta_unn[k])

        sIV.append((H ** 2.) * c[k] * Delta_nn[k])

        if A[k] != 0:

            P = c[k] * Delta_nn[k]
            p1 = (H ** 2.) * P
            p2 = (A_bar[k] ** 2.) * h[k] * Delta_unn[k]

            sV.append((A_bar[k] / A[k]) * (p1 + p2))
        else:
            sV.append(0)

        # ------------------------------------------------
        a = epsilon[k] * (c[k] ** 3.) * (h[k] ** 4.)

        sI_k.append(a * Delta_n[k] * ((h_bar[k] / h[k]) ** 0))

        sII_k.append(a * Delta_n[k] * ((h_bar[k] / h[k]) ** 1))

        sIII_k.append(a * Delta_n[k] * ((h_bar[k] / h[k]) ** 2))

        sIV_k.append(a * 0.)

        sV_k.append(a * Delta_n[k] * ((h_bar[k] / h[k]) ** 3))

        # ------------------------------------------------

    """ --------------------- Aspheric ----------------------------------------       
                      repeating for aspheric surf 
    
        
       Warren Smith - Modern optical engineering_ the design of optical systems
                       McGraw Hill (2008) pag 112
    ------------------------------------------------------------------------"""

    sI_a = []
    sII_a = []
    sIII_a = []
    sIV_a = []
    sV_a = []

    c = c_orig

    Asp = asp4 - ((asp2 / 4.0) * ((4.0 * asp2 ** 2) + (6.0 * asp2 * c) + (3.0 * c ** 2.0)))

    for k in range(0, len(c) - 1):
        a = 8.0 * Asp[k] * (h[k] ** 4.) * Delta_n[k]

        sI_a.append(a * ((h_bar[k] / h[k]) ** 0))

        sII_a.append(a * ((h_bar[k] / h[k]) ** 1))

        sIII_a.append(a * ((h_bar[k] / h[k]) ** 2))

        sIV_a.append(a * 0.)

        sV_a.append(a * ((h_bar[k] / h[k]) ** 3))

    """ ------------------------------------------------------------------
                               
                                   sums
    
    -----------------------------------------------------------------------"""

    print("Chief Ray Slope, Object Space : ", A_bar[0])
    print("Chief Ray Slope, Image Space : ", A_bar[-1])
    print("Marginal Ray Slope, Object Space : ", A[0])
    print("Marginal Ray Slope, Image Space : ", A[-1])
    print("Optical invariant: ", H)

    print("----------------------------------------------------------------")

    sI = np.asarray(sI)
    sII = np.asarray(sII)
    sIII = np.asarray(sIII)
    sIV = np.asarray(sIV)
    sV = np.asarray(sV)

    ###########################
    sI_k = np.asarray(sI_k)
    sII_k = np.asarray(sII_k)
    sIII_k = np.asarray(sIII_k)
    sIV_k = np.asarray(sIV_k)
    sV_k = np.asarray(sV_k)

    ###########################
    sI_a = np.asarray(sI_a)
    sII_a = np.asarray(sII_a)
    sIII_a = np.asarray(sIII_a)

    ###########################

    ### Seidel Aberration Coefficents
    si = sI - sI_k - sI_a
    sii = sII - sII_k - sII_a
    siii = sIII - sIII_k - sIII_a
    siv = sIV - sIV_k - sIV_a
    sv = sV - sV_k - sV_a
    
    SAC=[si, sii, siii, siv, sv]
    

    ### Seidel coefficients in waves
    Units = 1000.0
    W040 = Units * (1 / 8.0) * si / W
    W131 = Units * (1 / 2.0) * sii / W
    W222 = Units * (1 / 2.0) * siii / W
    W220 = Units * (1 / 4.0) * siv / W
    W311 = Units * (1 / 2.0) * sv / W
    
    SCW=[W040, W131, W222, W220, W311]
    
    ### Transverse Aberration Coefficents
    ## Spherical
    TSPH = si / (2.0 * u[-1] * n_1[-1])
    ## Coma
    TSCO = sii / (2.0 * u[-1] * n_1[-1])
    TTCO = 3.0 * sii / (2.0 * u[-1] * n_1[-1])
    ## Astigmatism
    TAST = -siii / (u[-1] * n_1[-1])
    ## Field Curvature
    TPFC = -siv / (2.0 * u[-1] * n_1[-1])
    TSFC = (siii + siv) / (2.0 * u[-1] * n_1[-1])
    TTFC = ((3.0 * siii) + siv) / (2.0 * u[-1] * n_1[-1])
    ## Distortion
    TDIS = -sv / (2.0 * u[-1] * n_1[-1])
    
    TAC=[TSPH, TSCO, TTCO, TAST, TPFC, TSFC, TTFC, TDIS]
    

    ### Longitudinal Aberration Coefficents
    ## Spherical
    LSPH = si / (2.0 * u[-1] * u[-1] * n_1[-1])
    ## Coma
    LSCO = sii / (2.0 * u[-1] * u[-1] * n_1[-1])
    LTCO = 3.0 * sii / (2.0 * u[-1] * u[-1] * n_1[-1])
    ## Astigmatism
    LAST = siii / (u[-1] * u[-1] * n_1[-1])
    ## Field Curvature
    LPFC = siv / (2.0 * u[-1] * u[-1] * n_1[-1])
    LSFC = (siii + siv) / (2.0 * u[-1] * u[-1] * n_1[-1])
    LTFC = ((3.0 * siii) + siv) / (2.0 * u[-1] * u[-1] * n_1[-1])
    ## Distortion
    LDIS = -sv / (2.0 * u[-1] * u[-1] * n_1[-1])
    LAC=[LSPH, LSCO, LTCO, LAST, LPFC, LSFC, LTFC, LDIS]

    SAC_N="Seidel Aberration Coefficents"
    SCW_N="Seidel coefficients in waves"
    TAC_N="Transverse Aberration Coefficents"
    LAC_N="Longitudinal Aberration Coefficents"
    
    AberNames=[SAC_N, SCW_N, TAC_N, LAC_N]
    return AberNames, SAC, SCW, TAC, LAC

