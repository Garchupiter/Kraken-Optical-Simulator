#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 13:30:17 2021

@author: joelherreravazquez
"""

import numpy as np
from .math_shapes_class import zernike__surf, zernike_expand, zernike_polynomials, zernike_math_notation


def RMS(SA, X, Y, Z, Zern_pol, z_pow):
    NZ = []
    SE = np.copy(SA)
    SE[0] = 0
    for i in range(np.shape(X)[0]):
        NZZ = Wavefront_Phase(X[i], Y[i], SE, Zern_pol, z_pow)
        NZ.append(NZZ)

    NZ = np.asarray(NZ)

    g = NZ - Z
    g = g * g
    g = np.mean(g)
    error = np.sqrt(g)

    gg = NZ
    gg = gg * gg
    gg = np.mean(gg)
    RMS = np.sqrt(gg)

    return RMS, error


def Zernike_Fitting(x1, y1, Z1, A, minimum=0.0001):

    NC = len(A)
    Zern_pol, z_pow = zernike_expand(NC)

    for i in range(0,2):
        Zi = System_Matrix_Zernikes(x1, y1, A, Zern_pol, z_pow, 0)
        ZT = Zi.T
        ZTZ = np.matmul(ZT, Zi)
        ZTZ_1 = np.linalg.inv(ZTZ)
        ZTZ_1_ZT = np.matmul(ZTZ_1, ZT)
        D = Z1
        D = np.asmatrix(D)
        D = D.T

        SA = np.zeros_like(A)
        NA = A.shape[0]
        MA = np.asarray(np.matmul(ZTZ_1_ZT, D))
        cont = 0

        ZZ = []
        for i in range(0, A.shape[0]):
            ZZ.append(zernike_math_notation(i, Zern_pol, z_pow))
        ZZ = np.asarray(ZZ)

        for i in range(0, NA):
            if A[i] != 0:
                SA[i] = MA[cont][0]

                cont = cont + 1
            else:
                SA[i] = 0.0

            # print(SA[i]," :   ",ZZ[[i][0]][:])


        A = np.abs(SA)
        Zeros = np.argwhere(A > minimum)
        AA = np.zeros_like(A)
        AA[Zeros] = 1
        A = AA



    WRMS, FITTINGERROR = RMS(SA, x1, y1, Z1, Zern_pol, z_pow)


    return SA, ZZ, WRMS#, FITTINGERROR






def Wavefront_Zernike_Phase(x, y, COEF):
    NC = len(COEF)
    Zern_pol, z_pow = zernike_expand(NC)

    tcoef = COEF.shape[0]
    p = np.sqrt((x * x) + (y * y))
    f = np.arctan2(x, y)  # angulo en direccion de las agujas del relog desde el eje x positivo

    ZFP = 0.0

    for i in range(0, tcoef):
        if COEF[i] != 0:
            ZFP = ZFP + COEF[i] * zernike_polynomials(i, p, f, Zern_pol, z_pow)
    return ZFP




def Wavefront_Phase(x, y, COEF, Zern_pol, z_pow):
    tcoef = COEF.shape[0]
    p = np.sqrt((x * x) + (y * y))
    f = np.arctan2(x, y)  # angulo en direccion de las agujas del relog desde el eje x positivo

    ZFP = 0.0

    for i in range(0, tcoef):
        if COEF[i] != 0:
            ZFP = ZFP + COEF[i] * zernike_polynomials(i, p, f, Zern_pol, z_pow)
    return ZFP


def Wf_XY_Components(x, y, N, Zern_pol, z_pow, fz):
    A = np.zeros(Zern_pol.shape[0])
    A[N] = 1.0
    z = Wavefront_Phase(x, y, A, Zern_pol, z_pow)

    return z


def System_Matrix_Zernikes(x, y, A, Zern_pol, z_pow, fz):
    NA = np.argwhere(A == 1)
    n_NA = NA.shape[0]

    Tp = x.shape[0]
    ZU = np.zeros((Tp, n_NA))

    cont = 0
    for h in range(0, Tp):
        for n in range(0, n_NA):
            F = Wf_XY_Components(x[h], y[h], NA[n], Zern_pol, z_pow, fz)

            ZU[cont, n] = F
        cont = cont + 1
    return ZU

#####################################################################################
