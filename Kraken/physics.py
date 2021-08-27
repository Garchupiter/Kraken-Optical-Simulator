# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 17:13:41 2021

@author: JOELHERRERAVAZQUEZ

"""
import numpy as np


############################################################################

def FresnelEnergy(vidrio, NP, NC, ImpVec, SurfNorm, ResVec, SETUP, Wave):
    """FresnelEnergy.

    :param vidrio:
    :param NP:
    :param NC:
    :param ImpVec:
    :param SurfNorm:
    :param ResVec:
    :param SETUP:
    :param Wave:
    """
    global Ts, Tp, Rs, Rp
    if vidrio != "MIRROR":
        Rp, Rs, Tp, Ts = fresnel_dielectric(NP, NC, ImpVec, SurfNorm, ResVec)

    if vidrio == "MIRROR":
        n_metal = np.interp(Wave, SETUP.W_alum, SETUP.N_alum)
        k_complex = np.interp(Wave, SETUP.W_alum, SETUP.K_alum)

        Rp, Rs, Tp, Ts = fresnel_metal(NP, n_metal, k_complex, ImpVec, SurfNorm)

    return Rp, Rs, Tp, Ts


##############################################################################

def fresnel_dielectric(NP, NC, LMN_Inc, LMN_nor_surf, LMN_result):
    """fresnel_dielectric.

    :param NP:
    :param NC:
    :param LMN_Inc:
    :param LMN_nor_surf:
    :param LMN_result:
    """
    CosTheta0 = np.abs(np.dot(LMN_Inc, LMN_nor_surf))
    CosTheta1 = np.abs(np.dot(LMN_result, LMN_nor_surf))
    n0 = NP
    n1 = NC

    #################### Enery ###################################
    # Chapter 5 - Handbook of Optics
    # Optical Society of America
    # Polarization
    # Section 5.6 pg.

    # rs, rp, ts,tp are the amplitudes reflection and refractions

    rs = (n0 * CosTheta0 - n1 * CosTheta1) / (n0 * CosTheta0 + n1 * CosTheta1)
    rp = (n1 * CosTheta0 - n0 * CosTheta1) / (n1 * CosTheta0 + n0 * CosTheta1)

    # No utilizado, no borrar- ts=2*n0*CosTheta0/(n0*CosTheta0+n1*CosTheta1)
    # No utilizado, no borrar - tp=2*n0*CosTheta0/(n1*CosTheta0+n0*CosTheta1)
    # Jones matrix for reflection [[-rp,0][0,rs]]
    # Jones matrix for transmition [[tp,0][0,ts]]
    # http://www.ub.edu/javaoptics/docs_applets/Doc_PolarEn.html

    # Intensity reflection R=r^2
    Rs = rs * rs
    Rp = rp * rp

    Tp = 1 - Rp
    Ts = 1 - Rs
    return Rp, Rs, Tp, Ts


##############################################################################
def fresnel_metal(NP, n_metal, k_complex, LMN_Inc, LMN_nor_surf):
    """fresnel_metal.

    :param NP:
    :param n_metal:
    :param k_complex:
    :param LMN_Inc:
    :param LMN_nor_surf:
    """

    n1 = NP
    n2 = np.complex(n_metal, k_complex)

    CosTheta0 = np.abs(np.dot(LMN_Inc, LMN_nor_surf))
    Thetai=np.arccos(CosTheta0)
    Thetat = np.arcsin(n1/n2*np.sin(Thetai))

    rs = (n1*np.cos(Thetai)-n2*np.cos(Thetat)) / (n1*np.cos(Thetai)+n2*np.cos(Thetat))
    rp = (n2*np.cos(Thetai)-n1*np.cos(Thetat)) / (n1*np.cos(Thetat)+n2*np.cos(Thetai))

    Rs = np.abs(rs)**2
    Rp = np.abs(rp)**2

    Tp = 1 - Rp
    Ts = 1 - Rs
    return Rp, Rs, Tp, Ts




######################################################################

def n_wave_dispersion(krakenSetup, GLSS, Wave):
    """n_wave_dispersion.

    :param krakenSetup:
    :param GLSS:
    :param Wave:
    """
    # allow direct refraction index asignment
    if GLSS[0:4] == "AIR_":
        n = float(GLSS[4::])
        Alpha = 0.0

    else:
        CAT = krakenSetup.CAT
        NAMES = krakenSetup.NAMES
        NM = krakenSetup.NM
        ED = krakenSetup.ED
        CD = krakenSetup.CD
        TD = krakenSetup.TD
        OD = krakenSetup.OD
        LD = krakenSetup.LD
        IT = krakenSetup.IT

        if GLSS == "MIRROR":
            # print("MIRROR")
            n = -1.  #
            Alpha = 0.0
        if GLSS == "AIR":
            # print("AIR")
            n = 1.  # n_air(Wave)
            Alpha = 0.0
        if GLSS != "MIRROR" and GLSS != "AIR":

            r = np.argwhere(NAMES == GLSS)[0][0]

            [Dispersion_Formula, MIL, Nd, Vd, Exclude_sub, Status] = NM[r]

            C = CD[r]  # Dispersion_Coefficents

            lam2 = Wave * Wave

            if Dispersion_Formula == 1:
                # print("Schott")
                a0, a1, a2, a3, a4, a5 = C[0], C[1], C[2], C[3], C[4], C[5]

                lam4 = lam2 * lam2
                lam6 = lam4 * lam2
                lam8 = lam6 * lam2

                n2 = a0 + (a1 * lam2) + (a2 / lam2) + (a3 / lam4) + (a4 / lam6) + (a5 / lam8)
                n = np.sqrt(n2)

            if Dispersion_Formula == 2:
                # print("Sellmeier 1")
                # print(GLSS)
                # print(C[0],C[1],C[2],C[3],C[4],C[5])

                K1, L1, K2, L2, K3, L3 = C[0], C[1], C[2], C[3], C[4], C[5]

                p1 = K1 * lam2 / (lam2 - L1)
                p2 = K2 * lam2 / (lam2 - L2)
                p3 = K3 * lam2 / (lam2 - L3)

                n = np.sqrt(p1 + p2 + p3 + 1.0)

            if Dispersion_Formula == 3:
                # print("Herzberger")
                # print(GLSS)
                # print(Dispersion_Formula,MIL,Nd,Vd,Exclude_sub,Status)

                CA = C[0]
                CB = C[1]
                CC = C[2]
                CD = C[3]
                CE = C[4]
                CF = C[5]

                CL = 1.0 / ((Wave * 2.) - 0.028)

                n = CA + (CB * CL) + (CC * CL ** 2.) + (CD * Wave ** 2.) + (CE * Wave ** 4.) + (CF * Wave ** 6.)

            if Dispersion_Formula == 4:
                print("Sellmeier 2")
            if Dispersion_Formula == 5:
                print("Conrady")
            if Dispersion_Formula == 6:
                print("Sellmeier 3")
            if Dispersion_Formula == 7:
                print("Hanbook of Optics 1")
            if Dispersion_Formula == 8:
                print("Hanbook of Optics 2")
            if Dispersion_Formula == 9:
                print("Sellmeier 4")
            if Dispersion_Formula == 10:
                print("Extended")
            if Dispersion_Formula == 11:
                print("Sellmeier 5")
            if Dispersion_Formula == 12:
                print("Extended 2")

            [Wa, Tr, Th] = IT[r]
            TR = np.interp(Wave, np.asarray(Wa), np.asarray(Tr))  # glass transmitance in specific wavelenght
            Alpha = -np.log(TR) / Th[0]  # Transmitance in the thickness
    # print("- - - - - - - - - - - - - - ")
    return n, Alpha


def ParaxCalc(N_Prec, SDT, SuTo, n, Gl):
    """ParaxCalc.

    :param N_Prec:
    :param SDT:
    :param SuTo:
    :param n:
    :param Gl:
    """
    # global RR, CurrN
    # Wave = WaveLength
    j = 0

    PrevN = N_Prec[j]

    S_Matrix = []
    N_Matrix = []

    CC = []

    # NN1.append(0)
    DD = []
    KK = []

    for j in range(0, n):

        sag_rad = SDT[j].Diameter / 400.0
        SG = SuTo.SurfaceShape(sag_rad, 0.0, j)

        # https://en.wikipedia.org/wiki/Sagitta_(optics)
        if SG != 0.0:
            R_p = (sag_rad * sag_rad) / (2.0 * SG)
        else:
            R_p = 999999999999.9

        # ------------------------------------------------

        if SDT[j].Solid_3d_stl != "None":
            R_p = 999999999999.9

        Glass = Gl[j]

        if Glass != "NULL":
            CurrN = N_Prec[j]

        if Glass == "NULL":
            CurrN = PrevN

        n1 = np.abs(PrevN)
        n2 = np.abs(CurrN)
        dd = SDT[j].Thickness + SDT[j].DespZ
        R = R_p

        if Glass == "MIRROR":
            n1 = -n1

        #     RR = np.matrix([[n1 / n2, (n1 - n2) / (n2 * R)], [0.0, 1]])

        #     # RR = np.matrix([[n1, -2 / R], [0.0, 1.0]])

        # if Glass != "MIRROR":

        #     RR = np.matrix([[n1 / n2, (n1 - n2) / (n2 * R)], [0.0, 1]])

        RR = np.matrix([[n1 / n2, (n1 - n2) / (n2 * R)], [0.0, 1]])

        CC.append(1 / R)

        # if NN1[-1]<0:
        #     NN1.append(-n1)

        # else:
        # NN2.append(CurrN)
        # NN1.append(PrevN)

        DD.append(dd)
        KK.append(0)

        TT = np.matrix([[1.0, 0.0], [dd, 1.0]])
        S_Matrix.append(RR)
        N_Matrix.append("R sup: " + str(j))

        S_Matrix.append(TT)
        N_Matrix.append("T sup: " + str(j) + " to " + str(j + 1))
        PrevN = CurrN

    SistemMatrix = np.matrix([[1.0, 0.0], [0.0, 1.0]])
    for Mi in S_Matrix[::-1]:
        SistemMatrix = np.dot(SistemMatrix, Mi)

    a = SistemMatrix[0][0, 0]
    b = SistemMatrix[0][0, 1]
    c = SistemMatrix[1][0, 0]
    d = SistemMatrix[1][0, 1]

    EFFL = -1 / b
    PPA = (1 - d) / b
    PPP = (a - 1) / b
    # print( "- - - - - - - - ")
    CC = np.asarray(CC)
    # print(CC)

    # print( "- - - - - - - - ")
    DD = np.asarray(DD)
    # print(DD)

    # print( "- - - - - - - - ")

    return SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, CC, N_Prec, DD
