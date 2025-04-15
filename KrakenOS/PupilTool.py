
import numpy as np
from .RayKeeper import *
from .Display import *
import scipy
from scipy import optimize
from scipy.optimize import fmin_cg
from scipy.optimize import fsolve
from .AstroAtmosphere import *
""" v2. """






import numpy as np
import matplotlib.pyplot as plt
import time
import math

from scipy.interpolate import interp1d
import scipy

import scipy
from scipy.optimize import minimize_scalar

##############################################################################

# Derivadas parciales para construir la matriz jacobiana
def jacobian(x, y, fun, h=1e-6):
    # Aproximación numérica de la jacobiana usando diferencias finitas
    f_x1 = fun(x + h, y)
    f_x2 = fun(x - h, y)
    f_y1 = fun(x, y + h)
    f_y2 = fun(x, y - h)

    j_11 = (f_x1[0] - f_x2[0]) / (2 * h)  # df1/dx
    j_12 = (f_y1[0] - f_y2[0]) / (2 * h)  # df1/dy
    j_21 = (f_x1[1] - f_x2[1]) / (2 * h)  # df2/dx
    j_22 = (f_y1[1] - f_y2[1]) / (2 * h)  # df2/dy

    return np.array([[j_11, j_12], [j_21, j_22]])

# Método de Newton-Raphson para minimizar rx y ry
def newton_raphson(x0, y0, fun, tol=1e-6, max_iter=20):
    x, y = x0, y0
    for i in range(max_iter):
        F = fun(x, y)
        J = jacobian(x, y, fun)

        # Verificar si la solución ha convergido
        if np.linalg.norm(F, ord=2) < tol:
            # print(f"Convergió en {i+1} iteraciones.")
            return x, y

        # Calcular el paso usando la inversa de la jacobiana
        try:
            delta = np.linalg.solve(J, F)  # J^-1 * F
        except np.linalg.LinAlgError:
            # print("La matriz jacobiana es singular.")
            return None

        # Actualizar las variables
        x, y = x - delta[0], y - delta[1]

    # print("No convergió dentro del número máximo de iteraciones.")
    x = 123456789.10
    y = 123456789.10
    return x, y




class RayStop:
    def __init__(self, System, Wave, StopSurf):
        """
        System : Sistema optico
        StopSurf : Superficie tomada como apertura del sistema
        Wave : Longitud de onda

        """

        self.System = System
        self.StopSurf = StopSurf
        self.Wave = Wave
        self.GlCorStop = np.asarray([0,0,0])
        self.System.IgnoreVignetting(0)
        self.IMP_LMN = [0,0,1]
        self.Ang_X = 0
        self.Ang_Y = 0
        self.rx = 0
        self.ry = 0
        self.preX = 0
        self.preY = 0


    def LMN(self, ix, iy):
        x0 = 0
        y0 = 0
        z0 = 0

        x1 = ix
        y1 = iy
        z1 = -10000.0

        r = np.sqrt((x1**2) + (y1**2) + (z1**2))

        L = (x1 - x0) / r
        M = (y1 - y0) / r
        N = (z1 - z0) / r

        return np.asarray([L, M, N])

    def LMN2(self, x0, y0, z0, x1, y1, z1):

        r = np.sqrt(((x1-x0)**2) + ((y1-y0)**2) + ((z1-z0)**2))

        L = (x1 - x0) / r
        M = (y1 - y0) / r
        N = (z1 - z0) / r

        return np.asarray([L, M, N])



    def CalculateField(self, V):

        ix =  V[0]
        iy =  V[1]
        dCos = self.LMN(ix, iy)
        pSource = self.GlCorStop

        self.System.RvTrace(pSource, dCos, self.Wave, self.StopSurf)
        xyz = self.System.XYZ

        [x0, y0, z0] = xyz[len(xyz)-1]
        [x1, y1, z1] = xyz[len(xyz)-2]

        Cx = x1 - x0
        Cy = y1 - y0
        Cz = z1 - z0

        AngX = np.rad2deg(np.arctan(Cx/Cz)) - self.Ang_X
        AngY = np.rad2deg(np.arctan(Cy/Cz)) - self.Ang_Y

        self.XYZ_I = [x0, y0, z0]
        self.XYZ_I2 = [x1, y1, z1]
        self.LMN_I = self.LMN2(x0, y0, z0, x1, y1, z1)

        return [AngX, AngY]

    def CalculateHeight(self, V):

        ix =  V[0]
        iy =  V[1]
        dCos = self.LMN(ix, iy)
        pSource = self.GlCorStop

        self.System.RvTrace(pSource, dCos, self.Wave, self.StopSurf)
        xyz = self.System.XYZ

        [x0, y0, z0] = xyz[len(xyz)-1]
        [x1, y1, z1] = xyz[len(xyz)-2]

        # Cx = x1 - x0
        # Cy = y1 - y0
        # Cz = z1 - z0

        # AngX = np.rad2deg(np.arctan(Cx/Cz)) - self.Ang_X
        # AngY = np.rad2deg(np.arctan(Cy/Cz)) - self.Ang_Y

        self.XYZ_I = [x0, y0, z0]
        self.XYZ_I2 = [x1, y1, z1]
        self.LMN_I = self.LMN2(x0, y0, z0, x1, y1, z1)

        return [x0- self.Ang_X, y0 - self.Ang_Y]


    def calculateXYZ(self, X, Y):

        dCos = self.INP_LMN
        pSource = [self.preX + X, self.preY + Y, 0]


        self.System.Trace(pSource, dCos, self.Wave)

        [x, y, z] = self.System.XYZ[-1]


        return np.array([x - self.Desp_X, y - self.Desp_Y])

    def calculateLMN(self, X, Y):

        [x0, y0, z0] = self.XYZ_I
        [x1, y1, z1] = self.XYZ_I2

        x1 = x1 + X
        y1 = y1 + Y


        r = np.sqrt(((x1 - x0)**2) + ((y1 - y0)**2) + ((z1 - z0)**2))

        L = (x1 - x0) / r
        M = (y1 - y0) / r
        N = (z1 - z0) / r

        dCos = [L,M, N]

        pSource = [self.preX, self.preY, 0]


        self.System.Trace(pSource, dCos, self.Wave)

        [x, y, z] = self.System.XYZ[-1]


        return np.array([x - self.Desp_X, y - self.Desp_Y])



def RMS_Pupil(r, SYSTEM, Surf, W, tet):
    """RMS_Pupil.

    Parameters
    ----------
    r :
        r
    SYSTEM :
        SYSTEM
    Surf :
        Surf
    W :
        W
    """
    SYSTEM.TargSurf(Surf)
    SYSTEM.IgnoreVignetting(0)
    SYSTEM.SurfFlat(Surf, 0)
    RP = raykeeper(SYSTEM)
    # tet = 0.01
    s_0 = [0.0, 0.0, 0.0]
    c_0 = [0.0, 0.0, 1.0]
    SYSTEM.Trace(s_0, c_0, W)
    RP.push()

    c_1 = [np.sin(np.deg2rad(tet)), 0.0, np.cos(np.deg2rad(tet))]
    s_0 = [r, 0.0, 0.0]
    SYSTEM.Trace(s_0, c_1, W)
    RP.push()

    c_2 = [np.sin(np.deg2rad((- tet))), 0.0, np.cos(np.deg2rad((- tet)))]
    s_0 = [(- r), 0.0, 0.0]
    SYSTEM.Trace(s_0, c_2, W)
    RP.push()

    c_3 = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
    s_0 = [0.0, r, 0.0]
    SYSTEM.Trace(s_0, c_3, W)
    RP.push()

    c_4 = [0.0, np.sin(np.deg2rad((- tet))), np.cos(np.deg2rad((- tet)))]
    s_0 = [0.0, (- r), 0.0]
    SYSTEM.Trace(s_0, c_4, W)
    RP.push()


    (X, Y, Z, L, M, N) = RP.pick(Surf)
    delta_Z = 0
    X = (((L / N) * delta_Z) + X)
    Y = (((M / N) * delta_Z) + Y)
    cenX = np.mean(X)
    cenY = np.mean(Y)
    x1 = (X - cenX)
    y1 = (Y - cenY)
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    SYSTEM.SurfFlat((- 1))
    SYSTEM.TargSurf((- 1))
    SYSTEM.Vignetting(0)
    RP.clean()
    return R_RMS

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

def FuncVectorCross(Z1, XYZa, LMNa, XYZb, LMNb, xy):
    """FuncVectorCross.

    Parameters
    ----------
    Z1 :
        Z1
    XYZa :
        XYZa
    LMNa :
        LMNa
    XYZb :
        XYZb
    LMNb :
        LMNb
    xy :
        xy
    """
    [La, Ma, Na] = LMNa
    [X0a, Y0a, Z0a] = XYZa
    [Lb, Mb, Nb] = LMNb
    [X0b, Y0b, Z0b] = XYZb
    if (xy == 0):
        V = ((((La / Na) * (Z1 - Z0a)) + X0a) - (((Lb / Nb) * (Z1 - Z0b)) + X0b))
    else:
        V = ((((Ma / Na) * (Z1 - Z0a)) + Y0a) - (((Mb / Nb) * (Z1 - Z0b)) + Y0b))
    return V

def DerVectCross(Z1, XYZa, LMNa, XYZb, LMNb, xy):
    """DerVectCross.

    Parameters
    ----------
    Z1 :
        Z1
    XYZa :
        XYZa
    LMNa :
        LMNa
    XYZb :
        XYZb
    LMNb :
        LMNb
    xy :
        xy
    """
    h = 1e-07
    f1 = FuncVectorCross((Z1 + h), XYZa, LMNa, XYZb, LMNb, xy)
    f2 = FuncVectorCross((Z1 - h), XYZa, LMNa, XYZb, LMNb, xy)
    der = ((f1 - f2) / (2 * h))
    return der

def SolveVectCross(XYZa, LMNa, XYZb, LMNb, xy):
    """SolveVectCross.

    Parameters
    ----------
    XYZa :
        XYZa
    LMNa :
        LMNa
    XYZb :
        XYZb
    LMNb :
        LMNb
    xy :
        xy
    """
    [La, Ma, Na] = LMNa
    [X0a, Y0a, Z0a] = XYZa
    [Lb, Mb, Nb] = LMNb
    [X0b, Y0b, Z0b] = XYZb
    Z1 = 1e-07
    cnt = 0
    while True:
        fun = FuncVectorCross(Z1, XYZa, LMNa, XYZb, LMNb, xy)
        derfun = DerVectCross(Z1, XYZa, LMNa, XYZb, LMNb, xy)
        Z2 = (Z1 - (fun / derfun))
        if (np.abs((Z1 - Z2)) < 1e-06):
            break
        else:
            Z1 = Z2
        if (cnt == 5):
            break
        cnt = (cnt + 1)
    X2 = (((La / Na) * (Z2 - Z0a)) + X0a)
    Y2 = (((Ma / Na) * (Z2 - Z0a)) + Y0a)
    return [X2, Y2, Z2]

def DerFpupil(SYSTEM, XY, H, Surf, tet, W, xy):
    """DerFpupil.

    Parameters
    ----------
    SYSTEM :
        SYSTEM
    XY :
        XY
    H :
        H
    Surf :
        Surf
    tet :
        tet
    W :
        W
    xy :
        xy
    """
    h = 1e-07
    der = (Fpupil(SYSTEM, (XY + h), H, Surf, tet, W, xy) - Fpupil(SYSTEM, (XY - h), H, Surf, tet, W, xy))
    DER = (der / (2.0 * h))
    return DER

def Fpupil(SYSTEM, XY, H, Surf, tet, W, xy):
    """Fpupil.

    Parameters
    ----------
    SYSTEM :
        SYSTEM
    XY :
        XY
    H :
        H
    Surf :
        Surf
    tet :
        tet
    W :
        W
    xy :
        xy
    """
    if (xy == 0):
        pSource_0 = [XY, 0.0, 0.0]
        dCos = [np.sin(np.deg2rad(tet)), 0.0, np.cos(np.deg2rad(tet))]
    else:
        pSource_0 = [0.0, XY, 0.0]
        dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
    SYSTEM.Trace(pSource_0, dCos, W)
    s = np.asarray(SYSTEM.SurfACE)
    a = np.squeeze(np.argwhere((s == Surf)))
    [X2, Y2, Z2] = SYSTEM.OST_XYZ[a]
    if (xy == 0):
        R = (X2 - H)
    else:
        R = (Y2 - H)
    return R

def SolveRayPupil(SYSTEM, H, tet, W, Surf, xy):
    """SolveRayPupil.

    Parameters
    ----------
    SYSTEM :
        SYSTEM
    H :
        H
    tet :
        tet
    W :
        W
    Surf :
        Surf
    xy :
        xy
    """
    XY0 = 1e-07
    cnt = 0
    while True:
        XY1 = (XY0 - (Fpupil(SYSTEM, XY0, H, Surf, tet, W, xy) / DerFpupil(SYSTEM, XY0, H, Surf, tet, W, xy)))
        if (np.abs((XY0 - XY1)) < 1e-06):
            break
        else:
            XY0 = XY1
        if (cnt == 5):
            break
        cnt = (cnt + 1)
    if (xy == 0):
        pSource_0 = [XY0, 0.0, 0.0]
        dCos = [np.sin(np.deg2rad(tet)), 0.0, np.cos(np.deg2rad(tet))]
    else:
        pSource_0 = [0.0, XY0, 0.0]
        dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]
    return (pSource_0, dCos)

class PupilCalc():
    """PupilCalc.
    """


    def __init__(self, system, Surf, W, ApTyp='STOP', AV=1.0):
        """__init__.

        Parameters
        ----------
        system :
            system
        Surf :
            Surf
        W :
            W
        ApTyp :
            ApTyp
        AV :
            AV
        """
        self.Surf = Surf
        self.W = W
        self.SYSTEM = system
        if (AV == 0):
            AV = 1.0
        self.FieldType = 'angle'
        self.ApertureType = ApTyp
        self.ApertureValue = AV
        self.x0 = 0
        self.y0 = 0
        self.z0 = 0
        self.L = 0.0
        self.M = 0.0
        self.N = 1.0
        self.Cordx = np.asarray(0)
        self.Cordy = np.asarray(0)
        self.Ptype = 'hexapolar'
        self.Samp = 6
        self.FieldX = 0.0
        self.FieldY = 0.0
        self.DirPupSal = [0.0, 0.0, 0.0]
        self.rad = 0
        self.theta = 0
        self.PupilInpFactor = 0.99
        self.menter = 1.0
        self.AtmosRef = 0
        self.T = 283.15
        self.P = 100500
        self.H = 0.0
        self.xc = 450
        self.lat = 50
        self.h = 0
        self.l1 = 5000.60169
        self.l2 = 0.50169
        self.z0 = 75.0
        self.tet= 0.01

        delta_Z = 0
        r = 1e-06
        cnt = 0
        while True:
            fun = RMS_Pupil(r, self.SYSTEM, Surf, self.W, self.tet)
            h = 1e-7
            f1 = RMS_Pupil((r + h), self.SYSTEM, Surf, self.W, self.tet)
            f2 = RMS_Pupil((r - h), self.SYSTEM, Surf, self.W, self.tet)
            der = ((f1 - f2) / (2 * h))
            r2 = (r - (fun / der))
            if (np.abs((r - r2)) < 1e-7):
                break
            else:
                r = r2
            if (cnt == 5):
                break
            cnt = (cnt + 1)
        self.SYSTEM.IgnoreVignetting(0)
        RP = raykeeper(self.SYSTEM)


        s_0 = [0.0, 0.0, 0.0]
        c_0 = [0.0, 0.0, 1.0]
        self.SYSTEM.Trace(s_0, c_0, self.W)

        RP.push()
        c_1 = [np.sin(np.deg2rad(self.tet)), 0.0, np.cos(np.deg2rad(self.tet))]
        s_0 = [r, 0.0, 0.0]
        self.SYSTEM.Trace(s_0, c_1, self.W)
        RP.push()

        c_2 = [np.sin(np.deg2rad((- self.tet))), 0.0, np.cos(np.deg2rad((- self.tet)))]
        s_0 = [(- r), 0.0, 0.0]
        self.SYSTEM.Trace(s_0, c_2, self.W)
        RP.push()

        c_3 = [0.0, np.sin(np.deg2rad(self.tet)), np.cos(np.deg2rad(self.tet))]
        s_0 = [0.0, r, 0.0]
        self.SYSTEM.Trace(s_0, c_3, self.W)
        RP.push()

        c_4 = [0.0, np.sin(np.deg2rad((- self.tet))), np.cos(np.deg2rad((- self.tet)))]
        s_0 = [0.0, (- r), 0.0]
        self.SYSTEM.Trace(s_0, c_4, self.W)
        RP.push()

# ----------------------------------------------------------------------------

        (X, Y, Z, L, M, N) = RP.pick(0)
        theta = 0

        for s in range(1, 5):
            theta = (theta + np.rad2deg(np.arccos((((L[0] * L[s]) + (M[0] * M[s])) + (N[0] * N[s])))))
        thetaIni = np.abs((theta / 4.0))

# ----------------------------------------------------------------------------

        (X, Y, Z, L, M, N) = RP.pick(Surf)
        theta = 0

        for s in range(1, 5):
            theta = (theta + np.rad2deg(np.arccos((((L[0] * L[s]) + (M[0] * M[s])) + (N[0] * N[s])))))
        thetaSurf = np.abs((theta / 4.0))

# ----------------------------------------------------------------------------

        (X, Y, Z, L, M, N) = RP.pick((- 1))
        theta = 0

        for s in range(1, 5):
            theta = (theta + np.rad2deg(np.arccos((((L[0] * L[s]) + (M[0] * M[s])) + (N[0] * N[s])))))
        thetaEnd = np.abs((theta / 4.0))

# ----------------------------------------------------------------------------



        M_ENTER_P = (thetaIni / thetaSurf)
        M_EXIT_P = (thetaIni / thetaEnd)
        STOP_DIAM = self.SYSTEM.SDT[Surf].Diameter

        if (self.ApertureType == 'EPD'):
            D_Input_Pup = self.ApertureValue
        if (self.ApertureType == 'STOP'):
            D_Input_Pup = (STOP_DIAM / M_ENTER_P)

        RadPupInp = (D_Input_Pup / 2.0)
        D_Exit_Pup = (STOP_DIAM * M_EXIT_P)
        RadPupOut = (D_Exit_Pup / 2.0)

        (Xs, Ys, Zs, Ls, Ms, Ns) = RP.pick(0)
        (Xe, Ye, Ze, Le, Me, Ne) = RP.pick((- 1))

        delta_Z = 0
        ZZ = (Ls, Ms, Ns, Xs, Ys)
        v0 = scipy.optimize.fsolve(R_RMS, delta_Z, args=ZZ)

        ZZ = (Le, Me, Ne, Xe, Ye)
        vf = scipy.optimize.fsolve(R_RMS, delta_Z, args=ZZ)

        PosPupInp = np.asarray([0, 0, v0[0]])
        OPS = np.asarray([Le[0], Me[0], Ne[0]])
        DirPupSal = OPS

        (X, Y, Z, L, M, N) = RP.pick((- 1))
        px = (((L[0] / N[0]) * vf) + X[0])
        py = (((M[0] / N[0]) * vf) + Y[0])
        pz = vf

        PosPupOut = np.asarray([px[0], py[0], (vf[0] + Ze[0])])
        PosPupOutFoc = np.asarray([px[0], py[0], pz[0]])

        self.SYSTEM.Vignetting(0)
        RP.clean()
        self.RadPupInp = RadPupInp
        self.PosPupInp = PosPupInp
        self.RadPupOut = RadPupOut
        self.PosPupOut = PosPupOut
        self.PosPupOutFoc = PosPupOutFoc
        self.DirPupSal = DirPupSal
        self.menter = M_ENTER_P
        Prx = self.SYSTEM.Parax(self.W)
        self.SistemMatrix, self.S_Matrix, self.N_Matrix, self.a, self.b, self.c, self.d, self.EFFL, self.PPA, self.PPP, self.CC, self.N_Prec, self.DD = Prx
        self.FocusAiryRadius=1.22*self.W*self.EFFL/(2.0*self.RadPupInp)






    def __patern_rect(self, x, y, kx, ky):
        """__patern_rect.

        Parameters
        ----------
        x :
            x
        y :
            y
        kx :
            kx
        ky :
            ky
        """
        for i in range((- (self.Samp * kx)), ((kx * self.Samp) + 1)):
            for j in range((- (self.Samp * ky)), ((ky * self.Samp) + 1)):
                x_0 = (i / self.Samp)
                y_0 = (j / self.Samp)
                r = np.sqrt(((x_0 * x_0) + (y_0 * y_0)))
                if (r <= 1.0):
                    x.append(x_0)
                    y.append(y_0)
        return (x, y)

    def Pattern(self):
        """Pattern.
        """
        self.Samp = int(self.Samp)
        x = []
        y = []
        if (self.Ptype == 'rtheta'):
            x.append((self.rad * np.cos(np.deg2rad(self.theta))))
            y.append((self.rad * np.sin(np.deg2rad(self.theta))))
        if (self.Ptype == 'chief'):
            x.append(0.0)
            y.append(0.0)
        if (self.Ptype == 'hexapolar'):
            x.append(0.0)
            y.append(0.0)
            for j in range(1, (self.Samp + 1)):
                m = (1.0 / self.Samp)
                r = (m * j)
                k = (j * 6)
                ang = (360.0 / k)
                for h in range(0, k):
                    x.append((r * np.cos(np.deg2rad((h * ang)))))
                    y.append((r * np.sin(np.deg2rad((h * ang)))))
        if (self.Ptype == 'square'):
            (x, y) = self.__patern_rect(x, y, 1, 1)
        if (self.Ptype == 'fanx'):
            (x, y) = self.__patern_rect(x, y, 1, 0)
        if (self.Ptype == 'fany'):
            (x, y) = self.__patern_rect(x, y, 0, 1)
        if (self.Ptype == 'fan'):
            (x, y) = self.__patern_rect(x, y, 1, 0)
            (x, y) = self.__patern_rect(x, y, 0, 1)
        if (self.Ptype == 'rand'):
            p = 1000000.0
            self.Sample = ((4 * self.Samp) * self.Samp)
            x_i = (np.random.randint((- p), p, self.Sample) / p)
            y_i = (np.random.randint((- p), p, self.Sample) / p)
            for i in range(0, self.Sample):
                x_0 = x_i[i]
                y_0 = y_i[i]
                r = np.sqrt(((x_0 * x_0) + (y_0 * y_0)))
                if (r < 1.0):
                    x.append(x_0)
                    y.append(y_0)
        self.Cordx = np.asarray(x)
        self.Cordy = np.asarray(y)

    def Pattern2Field(self):
        """Pattern2Field.
        """
        self.Pattern()
        x = (self.Cordx * self.RadPupInp)
        y = (self.Cordy * self.RadPupInp)

        (Px, Py, Pz) = self.PosPupInp
        if (self.FieldType == 'angle'):
            if (self.AtmosRef == 1):
                T = self.T
                P = self.P
                H = self.H
                xc = self.xc
                lat = self.lat
                h = self.h
                l1 = self.l1
                l2 = self.l2
                z0 = self.z0
                at = Observatory()
                n1 = at.n_tph(l=l1, T=T, p=P, RH=H, xc=xc)
                n2 = at.n_tph(l=l2, T=T, p=P, RH=H, xc=xc)
                rho = at.rho(p=P, T=T, RH=H, xc=xc)
                disp = dispersion(lat, h)
                disp.setReducedHeight(P, rho)
                f_x = (z0 + self.FieldX)
                f_y = self.FieldY
                Z0 = np.sqrt(((f_x ** 2) + (f_y ** 2)))
                atm_dispersion = disp.cassini(n1, n2, Z0)
                theta = np.arctan2(f_x, f_y)
                tx = (self.FieldX + (atm_dispersion * np.sin(theta)))
                ty = (self.FieldY + (atm_dispersion * np.cos(theta)))
                shiftX = (Pz * np.tan(np.deg2rad((- tx))))
                shiftY = (Pz * np.tan(np.deg2rad((- ty))))
            else:
                shiftX = (Pz * np.tan(np.deg2rad((- self.FieldX))))
                shiftY = (Pz * np.tan(np.deg2rad((- self.FieldY))))
            f_type = 1.0
        else:
            shiftX = (- self.FieldX)

            shiftY = (- self.FieldY)

            f_type = 0.0



        x0 = (np.copy((x * f_type)) + shiftX)
        y0 = (np.copy((y * f_type)) + shiftY)
        z = (np.ones_like(x) * Pz)
        z0 = 0.0*np.copy(x)
        X2 = ((x - x0) * (x - x0))
        Y2 = ((y - y0) * (y - y0))
        Z2 = ((z - z0) * (z - z0))
        S = np.sqrt(((X2 + Y2) + Z2))
        L = ((x - x0) / S)
        M = ((y - y0) / S)
        N = ((z - z0) / S)
        return (x0, y0, z0, L, M, N)




    def Pattern2FieldPlus(self):
        """Pattern2FieldPlus.
        """

        r = self.SYSTEM.SDT[self.Surf].Diameter/2.0

        X = self.Cordx * r
        Y = self.Cordy * r

        Z = np.zeros_like(X)
        V1 = np.ones_like(X)

        # print(X, Y, Z, " - . - . - . - . - . - . - . - . - . - . - . - ")


        StopPoint = np.array([X, Y, Z, V1])


        AA = self.SYSTEM.Pr3D.TRANS_1A[self.Surf].dot(StopPoint)


        XX = AA[0, :].tolist()
        YY = AA[1, :].tolist()
        ZZ = AA[2, :].tolist()
        VV = AA[3, :].tolist()

        RST = RayStop(self.SYSTEM, self.W, self.Surf)
        RST.Ang_X = self.FieldX
        RST.Ang_Y = self.FieldY

        H = 10000.0
        LimInf = [-H, -H]
        LimSup = [ H,  H]
        b=(LimInf, LimSup)

        Tx = 0
        Ty = 0

        # RST.GlCorStop = [-X[0][i], -Y[0][i], -Z[0][i]]

        StopPoint = [0, 0, 0, 1]
        AA = self.SYSTEM.Pr3D.TRANS_1A[self.Surf].dot(StopPoint)
        AA = AA.tolist()

        X = -AA[0][0]
        Y = -AA[0][1]
        Z = -AA[0][2]
        V = -AA[0][3]

        RST.GlCorStop = [X, Y, Z]
        R = scipy.optimize.least_squares(RST.CalculateField, [Tx,Ty],bounds = b,verbose = 0)

        xyz = RST.XYZ_I
        lmn = RST.LMN_I
        lmn = lmn.tolist()

        self.chief_xyz = xyz
        self.chief_lmn = lmn

        Tx = xyz[0]
        Ty = xyz[1]

        RST.System.IgnoreVignetting(0)
        RST.System.TargSurf(RST.StopSurf) # coloca como objetivo la superficie "Diafragma"

        RST.INP_LMN = lmn

        vx = []
        vy = []
        vz = []


        [ll,mm,nn] = lmn
        l = []
        m = []
        n = []



    ##############################################################################

        if (self.FieldType == 'angle'):

            for i in range(0, len(XX[0])):
                RST.preX = Tx # Coordenadas X del rayo en el origen
                RST.preY = Ty # Coordenadas Y del rayo en el origen

                RST.Desp_X = -XX[0][i] # Coordenada X deseada en el diafragma
                RST.Desp_Y = -YY[0][i] # Coordenada X deseada en el diafragma
                x0, y0 = 0, 0  # Punto inicial
                x, y = newton_raphson(x0, y0, RST.calculateXYZ)

                if x != 123456789.10:
                    vx.append(x + Tx)
                    vy.append(y + Ty)
                    vz.append(0)
                    l.append(ll)
                    m.append(mm)
                    n.append(nn)


    ##############################################################################

        if (self.FieldType == 'height'):

            for i in range(0, len(XX[0])):
                RST.preX = Tx # Coordenadas X del rayo en el origen
                RST.preY = Ty # Coordenadas Y del rayo en el origen

                RST.Desp_X = -XX[0][i] # Coordenada X deseada en el diafragma
                RST.Desp_Y = -YY[0][i] # Coordenada X deseada en el diafragma
                x0, y0 = 0, 0  # Punto inicial
                lx, my = newton_raphson(x0, y0, RST.calculateLMN)

                if lx != 123456789.10:
                    vx.append(Tx)
                    vy.append(Ty)
                    vz.append(0)

                    [x0, y0, z0] = RST.XYZ_I
                    [x1, y1, z1] = RST.XYZ_I2

                    x1 = x1 + lx
                    y1 = y1 + my

                    r = np.sqrt(((x1 - x0)**2) + ((y1 - y0)**2) + ((z1 - z0)**2))

                    L = (x1 - x0) / r
                    M = (y1 - y0) / r
                    N = (z1 - z0) / r

                    l.append(L)
                    m.append(M)
                    n.append(N)

        self.SYSTEM.TargSurfRest()
        return(vx, vy, vz, l, m, n)




