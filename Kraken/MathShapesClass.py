# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 17:37:31 2021

@author: JOELHERRERAVAZQUEZ
"""
import numpy as np
from scipy.interpolate import griddata

"""
Help
Zernike standard coefficents 1-37
#A[0]=0  #0                                                      Piston, Not Used
#A[1]=1 #1   4^(1/2) (p) * COS (A)    Tilt x,                   (about y axis)
#A[2]=1  #2   4^(1/2) (p) * SIN (A)    Tilt y,                   (about x axis)
#A[3]=1  #3   3^(1/2) (2p^2 - 1)                                  Power or Focus
#A[4]=1  #4   6^(1/2) (p^2) * SIN (2A)                          Astigmatism y, (45)
#A[5]=1  #5   6^(1/2) (p^2) * COS (2A)                          Astigmatism x, (0)
#A[6]=1  #6   8^(1/2) (3p^3 - 2p) * SIN (A)                      Coma y
#A[7]=1  #7   8^(1/2) (3p^3 - 2p) * COS (A)                      Coma x
#A[8]=1  #8   8^(1/2) (p^3) * SIN (3A)                          Trefoil y
#A[9]=1  #9   8^(1/2) (p^3) * COS (3A)                          Trefoil x
#A[10]=1 #10   5^(1/2) (6p^4 - 6p^2 + 1)                          Primary Spherical
#A[11]=1 #11  10^(1/2) (4p^4 - 3p^2) * COS (2A)                  Secondary Astigmatism x
#A[12]=1 #12 10^(1/2) (4p^4 - 3p^2) * SIN (2A)                  Secondary Astigmatism y
#A[13]=1 #13  10^(1/2) (p^4) * COS (4A)                          Tetrafoil x
#A[14]=1 #14  10^(1/2) (p^4) * SIN (4A)                          Tetrafoil y
#A[15]=1 #15  12^(1/2) (10p^5 - 12p^3 + 3p) * COS (A)          Secondary Coma x
#A[16]=1 #16  12^(1/2) (10p^5 - 12p^3 + 3p) * SIN (A)          Secondary Coma y
#A[17]=1 #17  12^(1/2) (5p^5 - 4p^3) * COS (3A)                  Secondary Trefoil x
#A[18]=1 #18 12^(1/2) (5p^5 - 4p^3) * SIN (3A)                  Secondary Trefoil y
#A[19]=1 #19  12^(1/2) (p^5) * COS (5A)                          Pentafoil x
#A[20]=1 #20  12^(1/2) (p^5) * SIN (5A)                          Pentafoil y
#A[21]=1 #21   7^(1/2) (20p^6 - 30p^4 + 12p^2 - 1)              Secondary Spherical
#A[22]=1 #22  14^(1/2) (15p^6 - 20p^4 + 6p^2) * SIN (2A)          Tertiary Astigmatism y
#A[23]=1 #23  14^(1/2) (15p^6 - 20p^4 + 6p^2) * COS (2A)          Tertiary Astigmatism x
#A[24]=1 #24  14^(1/2) (6p^6 - 5p^4) * SIN (4A)                  Secondary Tetrafoil y
#A[25]=1 #25  14^(1/2) (6p^6 - 5p^4) * COS (4A)                  Secondary Tetrafoil x
#A[26]=1 #26  14^(1/2) (p^6) * SIN (6A)
#A[27]=1 #27  14^(1/2) (p^6) * COS (6A)
#A[28]=1 #28  16^(1/2) (35p^7 - 60p^5 + 30p^3 - 4p) * SIN (A)  Tertiary Coma y
#A[29]=1 #29  16^(1/2) (35p^7 - 60p^5 + 30p^3 - 4p) * COS (A)  Tertiary Coma x
#A[30]=1 #30  16^(1/2) (21p^7 - 30p^5 + 10p^3) * SIN (3A)      Tertiary Trefoil y
#A[31]=1 #31  16^(1/2) (21p^7 - 30p^5 + 10p^3) * COS (3A)      Tertiary Trefoil x
#A[32]=1 #32  16^(1/2) (7p^7 - 6p^5) * SIN (5A)
#A[33]=1 #33  16^(1/2) (7p^7 - 6p^5) * COS (5A)
#A[34]=1 #34  16^(1/2) (p^7) * SIN (7A)
#A[35]=1 #35  16^(1/2) (p^7) * COS (7A)
#A[36]=1 #36   9^(1/2) (70p^8 - 140p^6 + 90p^4 - 20p^2 + 1)

__________________________
@author: JOELHERRERAVAZQUEZ
"""



class extra__surf:
    def __init__(self, C):
        """extra__surf."""

        """__init__.

        :param C:
        """
        self.COEF = C[1]
        self.user_surface=C[0]
    def calculate(self, x, y):
        """calculate.

        :param x:
        :param y:
        """
        Z = self.user_surface(x, y, self.COEF)
        return Z

##############################################################################
def even_asphere(x, y, E):
    """even_asphere.

    :param x:
    :param y:
    :param E:
    """
    r = np.sqrt((x * x) + (y * y))
    Z = np.zeros_like(x)
    for i in range(1, 9):
        if E[i - 1] != 0:
            Z = Z + E[i - 1] * np.power(r, 2 * i * 1.0)
            # Esto cae en un error para numeros muy grandes en potencias grandes
    return Z


##############################################################################
class aspheric__surf:
    def __init__(self, E):
        """aspheric__surf."""

        """__init__.

        :param E:
        """
        self.E = E

    def calculate(self, x, y):
        """calculate.

        :param x:
        :param y:
        """
        Z = even_asphere(x, y, self.E)
        return Z


##############################################################################
class conic__surf:
    def __init__(self, R_C, KON, C_RXY_RATIO):
        """conic__surf."""

        """__init__.

        :param R_C:
        :param KON:
        :param C_RXY_RATIO:
        """
        self.R_C = R_C
        self.KON = KON
        self.C_RXY_RATIO = C_RXY_RATIO

    def calculate(self, x, y):
        """calculate.

        :param x:
        :param y:
        """
        if self.R_C != 0.0:
            c = 1.0 / self.R_C
            s = np.sqrt((x * x) + ((y * self.C_RXY_RATIO) * (y * self.C_RXY_RATIO)))
            InRoot = 1 - (self.KON + 1.0) * c * c * s * s
            z = (c * s * s / (1.0 + np.sqrt(InRoot)))

        else:
            z = np.zeros_like(x)

        return z


###################################################################
class axicon__surf:
    def __init__(self, C_RXY_RATIO, AXC):
        """axicon__surf."""

        """__init__.

        :param C_RXY_RATIO:
        :param AXC:
        """

        self.C_RXY_RATIO = C_RXY_RATIO
        self.AXC = AXC

    def calculate(self, x, y):
        """calculate.

        :param x:
        :param y:
        """
        if self.AXC != 0:
            s = np.sqrt((x * x) + ((y * self.C_RXY_RATIO) * (y * self.C_RXY_RATIO)))
            z_axicon = s * np.tan(np.deg2rad(self.AXC))  # Axicon shape

        else:
            z_axicon = np.zeros_like(x)

        return z_axicon


###################################################################

class error_map__surf:
    def __init__(self, xValues, yValues, zValues, SPACE):
        """error_map__surf."""

        """__init__.

        :param xValues:
        :param yValues:
        :param zValues:
        :param SPACE:
        """

        TG = SPACE

        LGMX = np.max(yValues)
        LGMN = np.min(yValues)
        LG = LGMX - LGMN
        NPG = int(1 + LG / TG)

        s = 40
        VXX = np.arange((-LG / 2.0) - s * TG, (LG / 2.0) + s * TG, TG)
        VYY = np.arange((-LG / 2.0) - s * TG, (LG / 2.0) + s * TG, TG)

        NPG = VXX.shape[0]
        grid_x, grid_y = np.meshgrid(VXX, VYY)
        Z = np.zeros((NPG, NPG))

        err = 0.01
        cont = 1

        for h in range(0, NPG):
            for k in range(0, NPG):
                Ox = grid_x[h, k]
                Oy = grid_y[h, k]

                ARGWx = np.argwhere((xValues < (Ox + err)) & (xValues > (Ox - err)))
                ARGWy = np.argwhere((yValues < (Oy + err)) & (yValues > (Oy - err)))

                if (ARGWx.shape[0] > 0) and (ARGWy.shape[0] > 0):
                    IN = np.intersect1d(ARGWx, ARGWy)
                    if IN.shape[0] > 0:
                        varg = IN[0]
                        Z[h, k] = zValues[varg]

                        cont = cont + 1

        pnts = np.argwhere(Z != 0.0)
        values = Z[pnts[:, 0], pnts[:, 1]]
        Vx = VXX[pnts[:, 0]]
        Vy = VYY[pnts[:, 1]]
        points = np.vstack((Vx, Vy)).T
        aXX = np.arange((-LG / 2.0) - s * TG, (LG / 2.0) + s * TG, TG / 1.0)
        aYY = np.arange((-LG / 2.0) - s * TG, (LG / 2.0) + s * TG, TG / 1.0)
        grid_x, grid_y = np.meshgrid(aXX, aYY)

        grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')

        Z[0, ::] = grid_z0[0, ::]
        Z[NPG - 1, ::] = grid_z0[NPG - 1, ::]
        Z[::, 0] = grid_z0[::, 0]
        Z[::, NPG - 1] = grid_z0[::, NPG - 1]

        pnts = np.argwhere(Z != 0.0)

        Vx = VXX[pnts[:, 0]]
        Vy = VYY[pnts[:, 1]]

        self.values = Z[pnts[:, 0], pnts[:, 1]]
        self.points = np.vstack((Vx, Vy)).T

    # noinspection PyUnreachableCode
    def calculate(self, x, y):
        """calculate.

        :param x:
        :param y:
        """

        z = griddata(self.points, self.values, (x, y), method='cubic')
        return z
        23


###################################################################
class zernike__surf:

    def __init__(self, COEF, Z_POL, Z_POW, DMTR):
        """zernike__surf."""

        """__init__.

        :param COEF:
        :param Z_POL:
        :param Z_POW:
        :param DMTR:
        """
        self.COEF = COEF
        self.Z_POL = Z_POL
        self.Z_POW = Z_POW
        self.DMTR = DMTR

    def calculate(self, x, y):
        """calculate.

        :param x:
        :param y:
        """
        ZSP = np.zeros_like(x)
        for i in range(0, self.COEF.shape[0]):
            if self.COEF[i] != 0:
                p = (np.sqrt((x * x) + (y * y))) / (self.DMTR / 2.0)
                f = np.arctan2(x, y)  # angulo en direccion de las agujas del relog desde el eje x positivo
                ZSP = ZSP + self.COEF[i] * zernike_polynomials(i, p, f, self.Z_POL, self.Z_POW)
        return ZSP


def zernike_polynomials(term, ro, theta, Zern_pol, z_pow):
    """zernike_polynomials.

    :param term:
    :param ro:
    :param theta:
    :param Zern_pol:
    :param z_pow:
    """
    j, n, m, par, raiz = Zern_pol[term]
    ct = z_pow[term][0]
    pot = z_pow[term][1]
    NR = 0
    L = len(ct)
    for i in range(0, L):
        NR = ct[i] * np.power(ro, pot[i]) + NR
    if par == 1:
        S = raiz * NR
    if par == 3:
        S = raiz * NR * np.cos(m * theta)
    if par == 2:
        S = raiz * NR * np.sin(m * theta)
    return S


def z_parity(num):
    """z_parity.

    :param num:
    """
    nv = num / 2.0
    n = nv - int(nv)
    if n == 0:
        v = 2
    else:
        v = 3
    return v


def r_zern(m, n):
    """r_zern.

    :param m:
    :param n:
    """
    ls = int((n - m) / 2.0)
    cont = 0
    TCV = []
    pot = []
    a = []
    for s in range(0, int(ls) + 1):
        V1 = np.power(-1, s) * np.math.factorial(n - s)
        V2 = np.math.factorial(s) * np.math.factorial(((n + m) / 2.0) - s) * np.math.factorial(((n - m) / 2.0) - s)
        TC = V1 / V2
        potencia = (n - (2.0 * s))
        TCV.append(TC)
        pot.append(potencia)
        cont = cont + 1
    a.append(TCV)
    a.append(pot)
    return a


def zernike_expand(L):
    """zernike_expand.

    :param L:
    """
    cont = 0
    Z = np.array([0, 0, 0, 0, 0])
    n = L
    m = L
    for i in range(0, n):
        if cont >= L:
            break
        for j in range(0, m):
            v = (j - i) / 2.0
            if v - int(v) == 0:
                if i > j or i == j:
                    if j != 0:
                        v = z_parity(cont)
                        Z = np.vstack((Z, [cont, i, j, v, np.sqrt(2.0 * (i + 1.0))]))
                        if cont >= L:
                            break
                        cont = cont + 1
                        v = z_parity(cont)
                        Z = np.vstack((Z, [cont, i, j, v, np.sqrt(2.0 * (i + 1.0))]))
                        if cont >= L:
                            break
                        cont = cont + 1
                    if j == 0:
                        Z = np.vstack((Z, [cont, i, j, 1, np.sqrt(i + 1.0)]))
                        if cont >= L:
                            break
                        cont = cont + 1
    Z = np.delete(Z, 0, 0)
    E = []
    for i in range(0, Z.shape[0]):
        j, n, m, paR_c, raiz = Z[i]
        a = r_zern(m, n)
        E.append(a)
    Z = Z[:L]
    return Z, E


def zernike_math_notation(term, Zern_pol, z_pow):
    """zernike_math_notation.

    :param term:
    :param Zern_pol:
    :param z_pow:
    """
    j, n, m, par, raiz = Zern_pol[term]
    ZZZ = ["Piston",
           "Tilt x, (about y axis)",
           "Tilt y, (about x axis)",
           "Power or Focus",
           "Astigmatism y, (45deg)",
           "Astigmatism x, (0deg)",
           "Coma y",
           "Coma x",
           "Trefoil y",
           "Trefoil x",
           "Primary Spherical",
           "Secondary Astigmatism x",
           "Secondary Astigmatism y",
           "Tetrafoil x",
           "Tetrafoil y",
           "Secondary Coma x",
           "Secondary Coma y",
           "Secondary Trefoil x",
           "Secondary Trefoil y",
           "Pentafoil x",
           "Pentafoil y",
           "Secondary Spherical",
           "Tertiary Astigmatism y",
           "Tertiary Astigmatism x",
           "Secondary Tetrafoil y",
           "Secondary Tetrafoil x",
           "14^(1/2) (p^6) * SIN (6A)",
           "14^(1/2) (p^6) * COS (6A)",
           "Tertiary Coma y",
           "Tertiary Coma x",
           "Tertiary Trefoil y",
           "Tertiary Trefoil x",
           ]

    ct = z_pow[term][0]
    pot = z_pow[term][1]
    NR = " "
    L = len(ct)
    for i in range(0, L):
        if pot[i] == 0:
            NR = str(ct[i]) + "+" + NR
        else:
            NR = str(ct[i]) + "r^" + str(int(pot[i])) + "+" + NR
    i = len(NR)
    NR = NR[:-2]
    if m == 1:
        mm = ""
    else:
        mm = str(int(m))
    if par == 1:
        S = str(int(0.01 + raiz * raiz)) + "^(1/2)(" + NR + ")"
    if par == 3:
        S = str(int(0.01 + raiz * raiz)) + "^(1/2)(" + NR + ")" + "cos(" + mm + "T)"
    if par == 2:
        S = str(int(0.01 + raiz * raiz)) + "^(1/2)(" + NR + ")" + "sin(" + mm + "T)"
    if term < len(ZZZ):
        x = ZZZ[[term][0]][:]
    else:
        x = ""
    return S + "   " + x

############################################################################
