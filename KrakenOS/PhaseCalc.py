
import numpy as np
import KrakenOS as Kos
import scipy

def R_RMS(delta_Z, L, M, N, X, Y):
    """R_RMS.

    Parameters
    ----------
    delta_Z :
        delta_Z
    L :
        L
    M :
        M
    N :
        N
    X :
        X
    Y :
        Y
    """
    X = (((L / N) * delta_Z) + X)
    Y = (((M / N) * delta_Z) + Y)
    cenX = np.mean(X)
    cenY = np.mean(Y)
    x1 = (X - cenX)
    y1 = (Y - cenY)
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    return R_RMS

def Phase(PUPIL):
    """Phase.

    Parameters
    ----------
    PUPIL :
        PUPIL
    """
    SYSTEM = PUPIL.SYSTEM
    sup = PUPIL.Surf
    W = PUPIL.W
    ApType = PUPIL.ApertureType
    ApVal = PUPIL.ApertureValue
    configuracion_1 = SYSTEM.SETUP
    Samp = PUPIL.Samp
    Ptype = PUPIL.Ptype
    FieldY = PUPIL.FieldY
    FieldX = PUPIL.FieldX
    FieldType = PUPIL.FieldType
    pSource_0 = [0.0, 0.0, 0.0]
    dCos = [0.0, 0.0, 1.0]
    W = W
    RR = Kos.raykeeper(SYSTEM)
    RR2 = Kos.raykeeper(SYSTEM)
    Pup = Kos.PupilCalc(SYSTEM, sup, W, ApType, ApVal)
    Pup = Kos.PupilCalc(SYSTEM, sup, W, ApType, ApVal)
    Pup.Samp = Samp
    Pup.Ptype = Ptype
    Pup.FieldY = FieldY
    Pup.FieldX = FieldX
    Pup.FieldType = FieldType
    (x, y, z, L, M, N) = Pup.Pattern2Field()
    SYSTEM.IgnoreVignetting(0)
    for i in range(0, len(x)):
        pSource_0 = [x[i], y[i], z[i]]
        dCos = [L[i], M[i], N[i]]
        SYSTEM.Trace(pSource_0, dCos, W)
        RR.push()
        pSource_0 = (np.asarray([x[i], y[i], z[i]]) * 0.01)
        dCos = [L[i], M[i], N[i]]
        SYSTEM.Trace(pSource_0, dCos, W)
        RR2.push()
    (X_P, Y_P, Z_P, L_P, M_P, N_P) = RR2.pick((- 1))
    delta_Z = 0
    ZZ = (L_P, M_P, N_P, X_P, Y_P)
    v = scipy.optimize.fsolve(R_RMS, delta_Z, args=ZZ)
    X_P = (((L_P / N_P) * v) + X_P)
    Y_P = (((M_P / N_P) * v) + Y_P)
    CenX = np.mean(X_P)
    CenY = np.mean(Y_P)
    CenZ = (np.mean(Z_P) - v)
    parax_focus = np.asarray([CenX, CenY, CenZ[0]])
    (XS, YS, ZS, LS, MS, NS) = RR.pick((- 1))
    top = np.squeeze(RR.valid_TOP)
    op = np.squeeze(RR.valid_OP)
    op = op[:, (- 1)]
    (XPUP, YPUP, ZPUP, LPUP, MPUP, NPUP) = RR.pick(sup)
    RR.clean()
    RR2.clean()
    [Px, Py, Pz] = Pup.PosPupOut
    (X, Y, Z) = (CenX, CenY, CenZ[0])
    (x0, y0, z0) = (Px, Py, Pz)
    X2 = ((X - x0) * (X - x0))
    Y2 = ((Y - y0) * (Y - y0))
    Z2 = ((Z - z0) * (Z - z0))
    S = np.sqrt(((X2 + Y2) + Z2))
    L = ((X - x0) / S)
    M = ((Y - y0) / S)
    N = ((Z - z0) / S)
    OPSR = S
    Pup.Ptype = 'chief'
    (xc, yc, zc, Lc, Mc, Nc) = Pup.Pattern2Field()
    SYSTEM.IgnoreVignetting()
    pSource_0 = [xc[0], yc[0], zc[0]]
    dCos = [Lc[0], Mc[0], Nc[0]]
    SYSTEM.Trace(pSource_0, dCos, W)
    [Xc, Yc, Zc] = SYSTEM.XYZ[(- 1)]
    [Lc, Mc, Nc] = SYSTEM.LMN[(- 1)]
    XYZ0 = [Xc, Yc, Zc]
    LMN0 = [Lc, Mc, Nc]
    AP = S
    P_O = Kos.surf()
    P_O.Rc = 0
    P_O.Thickness = Pz
    P_O.Glass = 'AIR'
    P_O.Diameter = 1.0
    P_P = Kos.surf()
    SYSTEM.Parax(W)
    P_P.Rc = np.abs(AP)
    P_P.Thickness = 0
    P_P.Diameter = ((P_P.Rc * 0.99) * 2.0)
    P_P.Glass = 'MIRROR'
    P_P.DespX = Px
    P_P.DespY = Py
    P_P.TiltX = np.rad2deg(np.arcsin((- Mc)))
    P_P.TiltY = np.rad2deg(np.arcsin((Lc / np.cos(np.arcsin((- Mc))))))
    P_P.AxisMove = 0
    P_P.Order = 1
    P_I = Kos.surf()
    P_I.Diameter = (((P_P.Rc * 0.99) * 2.0) * 10.0)
    B = [P_O, P_P]
    PS = Kos.system(B, configuracion_1)
    RR = Kos.raykeeper(PS)
    PS.Trace(XYZ0, LMN0, W)
    OPSR = PS.TOP
    for i in range(0, len(XS)):
        pSource_0 = [XS[i], YS[i], ZS[i]]
        dCos = [LS[i], MS[i], NS[i]]
        PS.Trace(pSource_0, dCos, W)
        RR.push()
    VT = np.squeeze(RR.valid_TOP)
    VTOPSR = OPSR
    RRR = (((((VT - VTOPSR) * 1000.0) / W) * Pup.RadPupOut) / Pup.RadPupInp)
    RRRR = 1000.0
    P2V = (np.min(RRR) - np.max(RRR))
    return ((YPUP / 1000.0), (XPUP / 1000.0), RRR, P2V)

