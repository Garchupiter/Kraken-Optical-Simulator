# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Examp-2M-STL_ImageSlicer.py
"""

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


import numpy as np
import matplotlib.pyplot as plt
import scipy
import os

A1 = 1

if A1 == 0:
    # _________________________________________________________________#
    P_Obj = Kos.surf()
    P_Obj.Rc = 0
    P_Obj.Thickness = 1000 + 3.452200000000000E+003
    P_Obj.Glass = "AIR"
    P_Obj.Diameter = 1.059E+003 * 2.0
    # _________________________________________________________________#
    Thickness = 3.452200000000000E+003
    M1 = Kos.surf()
    M1.Rc = -9.638000000004009E+003
    M1.Thickness = -Thickness
    M1.k = -1.077310000000000E+000
    M1.Glass = "MIRROR"
    M1.Diameter = 1.059E+003 * 2.0
    M1.InDiameter = 250 * 2.0
    # _________________________________________________________________#
    M2 = Kos.surf()
    M2.Rc = -3.93E+003
    M2.Thickness = Thickness + 1.037525880125084E+003
    M2.k = -4.328100000000000E+000
    M2.Glass = "MIRROR"
    M2.Diameter = 3.365E+002 * 2.0
    M2.AxisMove = 0
    # _________________________________________________________________#
    P_Image_A = Kos.surf()
    P_Image_A.Diameter = 10.0
    P_Image_A.Glass = "AIR"
    P_Image_A.Thickness = 10
    P_Image_A.Name = "Image plane Tel"
    P_Image_A.DespZ = -100.0
    # _________________________________________________________________#
    A = [P_Obj, M1, M2, P_Image_A]

    # _________________________________________________________________#
    configuracion_1 = Kos.Setup()
    Telescopio = Kos.system(A, configuracion_1)
    Rayos = Kos.raykeeper(Telescopio)


    # _________________________________________________________________#

    # Gaussian
    def f(x):
        x = np.rad2deg(x)
        seing = 1.2 / 3600.0
        sigma = seing / 2.3548
        mean = 0
        standard_deviation = sigma
        y = scipy.stats.norm(mean, standard_deviation)
        res = y.pdf(x)
        return res


    Sun = Kos.SourceRnd()
    Sun.field = 4 * 1.2 / (2.0 * 3600.0)
    Sun.fun = f
    Sun.dim = 2100
    Sun.num = 100000
    L, M, N, X, Y, Z = Sun.rays()

    Xr = np.zeros_like(L)
    Yr = np.zeros_like(L)
    Zr = np.zeros_like(L)

    Lr = np.zeros_like(L)
    Mr = np.zeros_like(L)
    Nr = np.zeros_like(L)

    NM = np.zeros_like(L)

    con = 0
    con2 = 0
    W = 0.6

    for i in range(0, Sun.num):
        if con2 == 10:
            print(100. * i / Sun.num)
            con2 = 0

        pSource_0 = [X[i], Y[i], Z[i]]
        dCos = [L[i], M[i], N[i]]
        Telescopio.Trace(pSource_0, dCos, W)

        x, y, z = Telescopio.XYZ[-1]
        l, m, n = Telescopio.LMN[-1]
        Xr[con] = x
        Yr[con] = y
        Zr[con] = z
        Lr[con] = l
        Mr[con] = m
        Nr[con] = n

        if Telescopio.NAME[-1] == "Image plane Tel":
            NM[con] = i
        else:
            NM[con] = -1

        con = con + 1
        con2 = con2 + 1
        # Rayos.push()

    args = np.argwhere(NM != -1)

    X = Xr[args]
    Y = Yr[args]
    Z = Zr[args]

    L = Lr[args]
    M = Mr[args]
    N = Nr[args]
    W = W * np.ones_like(N)

    Rays = np.hstack((X, Y, Z, L, M, N, W))
    outfile = "savedRays.npy"
    np.save(outfile, Rays)









################################################################

else:
    P_Obj = Kos.surf()
    P_Obj.Rc = 0
    P_Obj.Thickness = 100. + 0.5
    P_Obj.Glass = "AIR"
    P_Obj.Diameter = 10

    currentDirectory = os.getcwd()
    direc = r"Jherrera-ImageSlicerBW-00.stl"
    P_ImageSlicer = Kos.surf()
    P_ImageSlicer.Diameter = 10.0
    P_ImageSlicer.Glass = "BK7"
    P_ImageSlicer.Name = "Image slicer"
    P_ImageSlicer.Solid_3d_stl = direc
    P_ImageSlicer.Thickness = 13
    P_ImageSlicer.TiltX = 180.0
    P_ImageSlicer.DespX = -0.55
    P_ImageSlicer.DespY = -0.03
    P_ImageSlicer.AxisMove = 0

    P_Ima = Kos.surf()
    P_Ima.Diameter = 10.0
    P_Ima.Glass = "AIR"
    P_Ima.Name = "Plano imagen"

    # _________________________________________________________________#

    A = [P_Obj, P_ImageSlicer, P_Ima]
    configuracion_1 = Kos.Setup()
    ImageSlicer = Kos.system(A, configuracion_1)
    Rayos = Kos.raykeeper(ImageSlicer)

    outfile = "savedRays.npy"
    R = np.load(outfile)
    print(np.shape(R))
    X, Y, Z, L, M, N, W = R[:, 0], R[:, 1], R[:, 2], R[:, 3], R[:, 4], R[:, 5], R[:, 6]

    nrays = 2000
    Xr = np.zeros(nrays)
    Yr = np.zeros(nrays)
    Zr = np.zeros(nrays)
    Lr = np.zeros(nrays)
    Mr = np.zeros(nrays)
    Nr = np.zeros(nrays)
    NM = np.zeros(nrays)

    con = 0
    con2 = 0

    for i in range(0, nrays):
        if con2 == 10:
            print(100. * i / nrays)
            con2 = 0

        pSource_0 = [X[i], Y[i], Z[i] * 0]
        dCos = [L[i], M[i], N[i]]
        ImageSlicer.NsTrace(pSource_0, dCos, W[i])

        x, y, z = ImageSlicer.XYZ[-1]
        l, m, n = ImageSlicer.LMN[-1]
        Xr[con] = x
        Yr[con] = y
        Zr[con] = z
        Lr[con] = l
        Mr[con] = m
        Nr[con] = n

        AA = ImageSlicer.SURFACE
        AA = np.asarray(AA)
        AW = np.argwhere(AA == 1)
        if ImageSlicer.NAME[-1] == "Plano imagen" and len(AW) < 10 and ImageSlicer.TT < 0.9:

            # and ImageSlicer.TT<0.9 and ImageSlicer.TT>0.4

            NM[con] = i
            Rayos.push()
        else:
            NM[con] = -1

        con = con + 1
        con2 = con2 + 1

    args = np.argwhere(NM != -1)

    X = Xr[args]
    Y = Yr[args]
    Z = Zr[args]

    L = Lr[args]
    M = Mr[args]
    N = Nr[args]
    W = W * np.ones_like(N)

    ###################
    plt.plot(X, Y, '.', c="r", markersize=1)

    # axis labeling
    plt.xlabel('x')
    plt.ylabel('y')

    # figure name
    plt.title('Dot Plot')
    plt.axis('square')
    plt.show()

    #             Rays.push()
    Kos.display3d(ImageSlicer, Rayos, 0)
