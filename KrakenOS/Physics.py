import numpy as np

def FresnelEnergy(vidrio, NP, NC, ImpVec, SurfNorm, ResVec, SETUP, Wave):
    """FresnelEnergy.

    Parameters
    ----------
    vidrio :
        vidrio
    NP :
        NP
    NC :
        NC
    ImpVec :
        ImpVec
    SurfNorm :
        SurfNorm
    ResVec :
        ResVec
    SETUP :
        SETUP
    Wave :
        Wave
    """
    # global Ts, Tp, Rs, Rp
    if (vidrio != 'MIRROR'):
        (Rp, Rs, Tp, Ts) = fresnel_dielectric(NP, NC, ImpVec, SurfNorm, ResVec)
    if (vidrio == 'MIRROR'):
        n_metal = np.interp(Wave, SETUP.W_alum, SETUP.N_alum)
        k_complex = np.interp(Wave, SETUP.W_alum, SETUP.K_alum)
        (Rp, Rs, Tp, Ts) = fresnel_metal(NP, n_metal, k_complex, ImpVec, SurfNorm)

    return (Rp, Rs, Tp, Ts)

def fresnel_dielectric(NP, NC, LMN_Inc, LMN_nor_surf, LMN_result):
    """fresnel_dielectric.

    Parameters
    ----------
    NP :
        NP
    NC :
        NC
    LMN_Inc :
        LMN_Inc
    LMN_nor_surf :
        LMN_nor_surf
    LMN_result :
        LMN_result
    """
    AA=(LMN_Inc[0]*LMN_nor_surf[0])+(LMN_Inc[1]*LMN_nor_surf[1])+(LMN_Inc[2]*LMN_nor_surf[2])
    CosTheta0 = np.abs(AA)
    BB=(LMN_result[0]*LMN_nor_surf[0])+(LMN_result[1]*LMN_nor_surf[1])+(LMN_result[2]*LMN_nor_surf[2])
    CosTheta1 = np.abs(BB)
    n0 = NP
    n1 = NC
    rs = (((n0 * CosTheta0) - (n1 * CosTheta1)) / ((n0 * CosTheta0) + (n1 * CosTheta1)))
    rp = (((n1 * CosTheta0) - (n0 * CosTheta1)) / ((n1 * CosTheta0) + (n0 * CosTheta1)))

    Rs = (rs ** 2.)
    Rp = (rp ** 2.)
    Tp = (1.0 - Rp)
    Ts = (1.0 - Rs)

    return (Rp, Rs, Tp, Ts)

def fresnel_metal(NP, n_metal, k_complex, LMN_Inc, LMN_nor_surf):
    """fresnel_metal.

    Parameters
    ----------
    NP :
        NP
    n_metal :
        n_metal
    k_complex :
        k_complex
    LMN_Inc :
        LMN_Inc
    LMN_nor_surf :
        LMN_nor_surf
    """
    n1 = NP
    n2 = np.complex(n_metal, k_complex)
    CosTheta0 = np.abs(np.dot(LMN_Inc, LMN_nor_surf))
    if CosTheta0> 1.0:

        CosTheta0 =1.0
    Thetai = np.arccos(CosTheta0)
    Thetat = np.arcsin(((n1 / n2) * np.sin(Thetai)))
    rs = (((n1 * np.cos(Thetai)) - (n2 * np.cos(Thetat))) / ((n1 * np.cos(Thetai)) + (n2 * np.cos(Thetat))))
    rp = (((n2 * np.cos(Thetai)) - (n1 * np.cos(Thetat))) / ((n1 * np.cos(Thetat)) + (n2 * np.cos(Thetai))))
    Rs = (np.abs(rs) ** 2)
    Rp = (np.abs(rp) ** 2)
    Tp = (1 - Rp)
    Ts = (1 - Rs)
    return (Rp, Rs, Tp, Ts)

def n_wave_dispersion(krakenSetup, GLSS, Wave):
    """n_wave_dispersion.

    Parameters
    ----------
    krakenSetup :
        krakenSetup
    GLSS :
        GLSS
    Wave :
        Wave
    """
    if (GLSS[0:4] != 'AIR_'):
        CAT = krakenSetup.CAT
        NAMES = krakenSetup.NAMES
        NM = krakenSetup.NM
        ED = krakenSetup.ED
        CD = krakenSetup.CD
        TD = krakenSetup.TD
        OD = krakenSetup.OD
        LD = krakenSetup.LD
        IT = krakenSetup.IT
        if (GLSS == 'MIRROR'):
            n = (- 1.0)
            Alpha = 0.0

        if (GLSS == 'AIR'):
            n = 1.0
            Alpha = 0.0
        if ((GLSS != 'MIRROR') and (GLSS != 'AIR')):
            # print(GLSS)
            r = np.argwhere((NAMES == GLSS))[0][0]
            [Dispersion_Formula, MIL, Nd, Vd, Exclude_sub, Status] = NM[r]
            C = CD[r]
            lam2 = Wave **2.

            if (Dispersion_Formula == 1):
                """ Schott """
                (a0, a1, a2, a3, a4, a5) = (C[0], C[1], C[2], C[3], C[4], C[5])
                lam4 = Wave **4.
                lam6 = Wave **6.
                lam8 = Wave **8.
                n2 = (((((a0 + (a1 * lam2)) + (a2 / lam2)) + (a3 / lam4)) + (a4 / lam6)) + (a5 / lam8))
                n = np.sqrt(n2)

            if (Dispersion_Formula == 2):
                """ Sellmeier1 """
                (K1, L1, K2, L2, K3, L3) = (C[0], C[1], C[2], C[3], C[4], C[5])
                p1 = ((K1 * lam2) / (lam2 - L1))
                p2 = ((K2 * lam2) / (lam2 - L2))
                p3 = ((K3 * lam2) / (lam2 - L3))
                n = np.sqrt((((p1 + p2) + p3) + 1.0))

            if (Dispersion_Formula == 3):
                """Herzberger"""
                CA = C[0]
                CB = C[1]
                CC = C[2]
                CD = C[3]
                CE = C[4]
                CF = C[5]
                CL = (1.0 / ((Wave * 2.0) - 0.028))
                n = (((((CA + (CB * CL)) + (CC * (CL ** 2.0))) + (CD * (Wave ** 2.0))) + (CE * (Wave ** 4.0))) + (CF * (Wave ** 6.0)))
            if (Dispersion_Formula == 4):
                """Sellmeier 2"""
                GLSN = C[0] + (C[1] * Wave**2 / (Wave**2 - (C[2])**2)) + (C[3] * Wave**2 / (Wave**2 - (C[4])**2))
                n = np.sqrt(GLSN + 1.0)

            if (Dispersion_Formula == 5):
                """Conrady"""
                n = C[0] + (C[1] / Wave) + (C[2] / Wave**3.5)

            if (Dispersion_Formula == 6):
                """Sellmeier 3"""
                GLSN = (C[0] * Wave**2 / (Wave**2 - C[1])) + (C[2] * Wave**2 / (Wave**2 - C[3])) + (C[4] * Wave**2 / (Wave**2 - C[5])) + (C[6] * Wave**2 / (Wave**2 - C[7]))
                n = np.sqrt(GLSN + 1.0)

            if (Dispersion_Formula == 7):
                """Hanbook of Optics 1"""
                GLSN = C[0] + (C[1] / (Wave**2 - C[2])) - (C[3] * Wave**2)
                n = np.sqrt(GLSN)

            if (Dispersion_Formula == 8):
                """Hanbook of Optics 2"""
                GLSN = C[0] + (C[1] * Wave**2 / (Wave**2 - C[2])) - (C[3] * Wave**2)
                n = np.sqrt(GLSN)

            if (Dispersion_Formula == 9):
                """Sellmeier 4"""
                GLSN = C[0] + (C[1] * Wave**2 / (Wave**2 - C[2])) + (C[3] * Wave**2 / (Wave**2 - C[4]))
                n = np.sqrt(GLSN)

            if (Dispersion_Formula == 10):
                """Extended"""
                GLSN = C[0] + (C[1] * Wave**2) + (C[2] * Wave**-2) + (C[3] * Wave**-4) + (C[4] * Wave**-6) + \
                              (C[5] * Wave**-8) + (C[6] * Wave**-10) + (C[7] * Wave**-12)
                n = np.sqrt(GLSN)

            if (Dispersion_Formula == 11):
                """Sellmeier 5"""
                GLSN = (C[0] * Wave**2 / (Wave**2 - C[1])) + (C[2] * Wave**2 / (Wave**2 - C[3])) + (C[4] * Wave**2 / (Wave**2 - C[5])) + (C[6] * Wave**2 / (Wave**2 - C[7])) + (C[8] * Wave**2 / (Wave**2 - C[9]))
                n = np.sqrt(GLSN + 1.0)


            if (Dispersion_Formula == 12):
                """Extended 2"""
                GLSN = C[0] + (C[1] * Wave**2) + (C[2] * Wave**-2) + (C[3] * Wave**-4) + (C[4] * Wave**-6) + (C[5] * Wave**-8) + (C[6] * Wave**4) + (C[7] * Wave**6)
                n = np.sqrt(GLSN)


            if (Dispersion_Formula == 13):
                """HIKARI Catalog"""
                GLSN0 = C[0]
                GLSN1 =(C[1] * Wave**2)
                GLSN2 =(C[2] * Wave**4)
                GLSN3 =(C[3] * Wave**-2)
                GLSN4 =(C[4] * Wave**-4)
                GLSN5 =(C[5] * Wave**-6)
                GLSN6 =(C[6] * Wave**-8)
                GLSN7 =(C[7] * Wave**-10)
                GLSN8 =(C[8] * Wave**-12)
                GLSN = GLSN0 + GLSN1 + GLSN2 + GLSN3 + GLSN4 + GLSN5 + GLSN6 + GLSN7 + GLSN8

                """ n(W)2 = A0 + A1*W**2 + A2*W**4 + A3*W**−2 + A4*W**−4 + A5*W**−6 + A6*W**−8 + A7*W**−10 + A8*W**−12"""

                n = np.sqrt(GLSN)

            [Wa, Tr, Th] = IT[r]
            TR = np.interp(Wave, np.asarray(Wa), np.asarray(Tr))
            Alpha = ((- np.log(TR)) / Th[0])

    else:
        n = float(GLSS[4:])
        Alpha = 0.0

    return (n, Alpha)

def ParaxCalc(N_Prec, SDT, SuTo, n, Gl):
    """ParaxCalc.

    Parameters
    ----------
    N_Prec :
        N_Prec
    SDT :
        SDT
    SuTo :
        SuTo
    n :
        n
    Gl :
        Gl
    """
    j = 0
    PrevN = N_Prec[j]
    S_Matrix = []
    N_Matrix = []
    CC = []
    DD = []
    KK = []
    for j in range(0, n+0):
        sag_rad = (SDT[j].Diameter / 400.0)
        SG = SuTo.SurfaceShape(sag_rad, 0.0, j)
        if (SG != 0.0):
            R_p = ((sag_rad * sag_rad) / (2.0 * SG))
        else:
            R_p = 9999999999999.9
        if (SDT[j].Solid_3d_stl != 'None'):
            R_p = 9999999999999.9
        Glass = Gl[j]
        if (Glass != 'NULL'):
            CurrN = N_Prec[j]
        if (Glass == 'NULL'):
            CurrN = PrevN
        n1 = np.abs(PrevN)
        n2 = np.abs(CurrN)
        dd = (SDT[j].Thickness + SDT[j].DespZ)
        R = R_p
        if (Glass == 'MIRROR'):
            n1 = (- n1)

        if (SDT[j].Thin_Lens == 0):
            RR = np.matrix([[(n1 / n2), ((n1 - n2) / (n2 * R))], [0.0, 1]])
        else:
            RR = np.matrix([[1.0, -1.0 / SDT[j].Thin_Lens], [0.0, 1]])

        CC.append((1 / R))
        DD.append(dd)
        KK.append(0)
        TT = np.matrix([[1.0, 0.0], [dd, 1.0]])
        S_Matrix.append(RR)
        N_Matrix.append(('R surf ' + str(j)))
        S_Matrix.append(TT)
        N_Matrix.append(((('T surf ' + str(j)) + ' to ') + str((j + 1))))
        PrevN = CurrN
    SistemMatrix = np.matrix([[1.0, 0.0], [0.0, 1.0]])
    for Mi in S_Matrix[::(- 1)]:
        SistemMatrix = np.dot(SistemMatrix, Mi)

    """ Hecht E. - Optics-AW (2002) """
    a = SistemMatrix[0][(0, 0)]
    b = SistemMatrix[0][(0, 1)]
    c = SistemMatrix[1][(0, 0)]
    d = SistemMatrix[1][(0, 1)]
    EFFL = ((- 1) / b)
    PPA = ((1 - a) / -b)
    PPP = ((d - 1) / -b)
    CC = np.asarray(CC)
    DD = np.asarray(DD)
    return (SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, CC, N_Prec, DD)











