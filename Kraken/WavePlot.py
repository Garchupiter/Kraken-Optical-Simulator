#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 20:13:17 2021

@author: joelherreravazquez
"""

import numpy as np
import matplotlib.pyplot as plt
import Kraken as kn


def ZernikeDataImage2Plot(datos, Type="interferogram"):
    """ZernikeDataImage2Plot.

    :param datos:
    :param Type:
    """

    if Type == "interferogram":
        datos = np.sin(2 * np.pi * datos)
    if Type == "phase":
        datos = datos

    img = datos / np.max(datos)
    img = img * 255
    img = img.astype(int)
    img = np.asarray(img)
    plt.figure(987)
    plt.imshow(img, cmap='gray')
    plt.show()
    return 0


def WavefrontData2Image(z_coeff, res=323):
    """WavefrontData2Image.

    :param z_coeff:
    :param res:
    """
    TamImag = int(res)
    r = TamImag / 2.0
    ARRAY_ZERNIKE = np.zeros((TamImag, TamImag))

    H = []
    K = []
    X = []
    Y = []
    for h in range(0, TamImag):
        for k in range(0, TamImag):
            x = (h - r) / r
            y = (k - r) / r

            RP = np.sqrt((x * x) + (y * y))

            if RP <= 1:
                H.append(h)
                K.append(k)
                X.append(x)
                Y.append(y)

    H = np.asarray(H)
    K = np.asarray(K)
    X = np.asarray(X)
    Y = np.asarray(Y)

    Z = kn.Wavefront_Zernike_Phase(X, Y, z_coeff)
    ARRAY_ZERNIKE[H, K] = Z

    return ARRAY_ZERNIKE
