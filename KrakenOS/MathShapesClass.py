import numpy as np
from scipy.interpolate import griddata

class extra__surf():
    """extra__surf.
    """


    def __init__(self, C):
        """__init__.

        Parameters
        ----------
        C :
            C
        """
        self.COEF = C[1]
        self.user_surface = C[0]

    def calculate(self, x, y):
        """calculate.

        Parameters
        ----------
        x :
            x
        y :
            y
        """
        Z = self.user_surface(x, y, self.COEF)
        return Z

def even_asphere(x, y, E):
    """even_asphere.

    Parameters
    ----------
    x :
        x
    y :
        y
    E :
        E
    """
    r = np.sqrt(((x * x) + (y * y)))
    Z = 0.0 * np.zeros_like(x)
    for i in range(1, 9):
        if (E[(i - 1)] != 0):
            Z = (Z + (E[(i - 1)] * np.power(r, ((2.0 * i*1.0) * 1.0))))
    return Z

class aspheric__surf():
    """aspheric__surf.
    """
    def __init__(self, E):
        """__init__.

        Parameters
        ----------
        E :
            E
        """
        self.E = E

    def calculate(self, x, y):
        """calculate.

        Parameters
        ----------
        x :
            x
        y :
            y
        """
        Z = even_asphere(x, y, self.E)
        return Z

class conic__surf(object):
    """conic__surf.
    """
    def __init__(self, R_C, KON, C_RXY_RATIO):
        """__init__.

        Parameters
        ----------
        R_C :
            R_C
        KON :
            KON
        C_RXY_RATIO :
            C_RXY_RATIO
        """
        self.R_C = R_C
        self.KON = KON
        self.C_RXY_RATIO = C_RXY_RATIO

    def calculate(self, x, y):
        """calculate.

        Parameters
        ----------
        x :
            x
        y :
            y
        """

        z = CalculateCon(x, y, self.R_C , self.C_RXY_RATIO, self.KON)

        return z

def CalculateCon(x, y, R_C , C_RXY_RATIO, KON):
    """calculate.

    Parameters
    ----------
    x :
        x
    y :
        y
    """
    if (R_C != 0.0):

        s = np.sqrt(((x * x) + ((y * C_RXY_RATIO) * (y * C_RXY_RATIO))))

        c = (1.0 / R_C)
        InRoot = (1 - (((((KON + 1.0) * c) * c) * s) * s))

        InRoot = np.abs(InRoot)
        z = (((c * s) * s) / (1.0 + np.sqrt(InRoot)))
    else:
        z = 0.0 * np.zeros_like(x)

    return z

class axicon__surf():
    """axicon__surf.
    """

    def __init__(self, C_RXY_RATIO, AXC):
        """__init__.

        Parameters
        ----------
        C_RXY_RATIO :
            C_RXY_RATIO
        AXC :
            AXC
        """
        self.C_RXY_RATIO = C_RXY_RATIO
        self.AXC = AXC

    def calculate(self, x, y):
        """calculate.

        Parameters
        ----------
        x :
            x
        y :
            y
        """

        z_axicon = CalculateAxic( x, y, self.C_RXY_RATIO, self.AXC)
        return z_axicon

def CalculateAxic( x, y, C_RXY_RATIO, AXC):
    """calculate.

    Parameters
    ----------
    x :
        x
    y :
        y
    """
    if (AXC != 0.0):
        s = np.sqrt(((x * x) + ((y * C_RXY_RATIO) * (y * C_RXY_RATIO))))
        z_axicon = (s * np.tan(np.deg2rad(AXC)))
    else:
        z_axicon = 0.0 * np.zeros_like(x)
    return z_axicon

class error_map__surf():
    """error_map__surf.
    """

    def __init__(self, xValues, yValues, zValues, SPACE):
        """__init__.

        Parameters
        ----------
        xValues :
            xValues
        yValues :
            yValues
        zValues :
            zValues
        SPACE :
            SPACE
        """
        TG = SPACE
        LGMX = np.max(yValues)
        LGMN = np.min(yValues)
        LG = (LGMX - LGMN)
        NPG = int((1 + (LG / TG)))
        s = 40
        VXX = np.arange((((- LG) / 2.0) - (s * TG)), ((LG / 2.0) + (s * TG)), TG)
        VYY = np.arange((((- LG) / 2.0) - (s * TG)), ((LG / 2.0) + (s * TG)), TG)
        NPG = VXX.shape[0]
        (grid_x, grid_y) = np.meshgrid(VXX, VYY)
        Z = np.zeros((NPG, NPG))
        err = 0.01
        cont = 1
        for h in range(0, NPG):
            for k in range(0, NPG):
                Ox = grid_x[(h, k)]
                Oy = grid_y[(h, k)]
                ARGWx = np.argwhere(((xValues < (Ox + err)) & (xValues > (Ox - err))))
                ARGWy = np.argwhere(((yValues < (Oy + err)) & (yValues > (Oy - err))))
                if ((ARGWx.shape[0] > 0) and (ARGWy.shape[0] > 0)):
                    IN = np.intersect1d(ARGWx, ARGWy)
                    if (IN.shape[0] > 0):
                        varg = IN[0]
                        Z[(h, k)] = zValues[varg]
                        cont = (cont + 1)
        pnts = np.argwhere((Z != 0.0))
        values = Z[(pnts[:, 0], pnts[:, 1])]
        Vx = VXX[pnts[:, 0]]
        Vy = VYY[pnts[:, 1]]
        points = np.vstack((Vx, Vy)).T
        aXX = np.arange((((- LG) / 2.0) - (s * TG)), ((LG / 2.0) + (s * TG)), (TG / 1.0))
        aYY = np.arange((((- LG) / 2.0) - (s * TG)), ((LG / 2.0) + (s * TG)), (TG / 1.0))
        (grid_x, grid_y) = np.meshgrid(aXX, aYY)
        grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
        Z[0, :] = grid_z0[0, :]
        Z[(NPG - 1), :] = grid_z0[(NPG - 1), :]
        Z[:, 0] = grid_z0[:, 0]
        Z[:, (NPG - 1)] = grid_z0[:, (NPG - 1)]
        pnts = np.argwhere((Z != 0.0))
        Vx = VXX[pnts[:, 0]]
        Vy = VYY[pnts[:, 1]]
        self.values = Z[(pnts[:, 0], pnts[:, 1])]
        self.points = np.vstack((Vx, Vy)).T

    def calculate(self, x, y):
        """calculate.

        Parameters
        ----------
        x :
            x
        y :
            y
        """
        z = griddata(self.points, self.values, (x, y), method='cubic')
        return z

class zernike__surf():
    """zernike__surf.
    """

    def __init__(self, COEF, Z_POL, Z_POW, DMTR):
        """__init__.

        Parameters
        ----------
        COEF :
            COEF
        Z_POL :
            Z_POL
        Z_POW :
            Z_POW
        DMTR :
            DMTR
        """
        self.COEF = COEF
        self.Z_POL = Z_POL
        self.Z_POW = Z_POW
        self.DMTR = DMTR

    def calculate(self, x, y):
        """calculate.

        Parameters
        ----------
        x :
            x
        y :
            y
        """

        ZSP = CalculateZern(x, y, self.Z_POL, self.Z_POW, self.COEF, self.DMTR)
        return ZSP

def CalculateZern( x, y, Z_POL, Z_POW, COEF, DMTR):
    """calculate.

    Parameters
    ----------
    x :
        x
    y :
        y
    """
    ZSP = 0.0*np.zeros_like(x)
    for i in range(0, COEF.shape[0]):
        if (COEF[i] != 0):
            p = (np.sqrt(((x * x) + (y * y))) / (DMTR / 2.0))
            f = np.arctan2(x, y)
            ZSP = (ZSP + (COEF[i] * zernike_polynomials(i, p, f, Z_POL, Z_POW)))
    return ZSP

def zernike_polynomials(term, ro, theta, Zern_pol, z_pow):
    """zernike_polynomials.

    Parameters
    ----------
    term :
        term
    ro :
        ro
    theta :
        theta
    Zern_pol :
        Zern_pol
    z_pow :
        z_pow
    """
    (j, n, m, par, raiz) = Zern_pol[term]
    ct = z_pow[term][0]
    pot = z_pow[term][1]
    NR = 0
    L = len(ct)
    for i in range(0, L):
        NR = ((ct[i] * np.power(ro, pot[i])) + NR)
    if (par == 1):
        S = (raiz * NR)
    if (par == 3):
        S = ((raiz * NR) * np.cos((m * theta)))
    if (par == 2):
        S = ((raiz * NR) * np.sin((m * theta)))
    return S

def z_parity(num):
    """z_parity.

    Parameters
    ----------
    num :
        num
    """
    nv = (num / 2.0)
    n = (nv - int(nv))
    if (n == 0):
        v = 2
    else:
        v = 3
    return v

def r_zern(m, n):
    """r_zern.

    Parameters
    ----------
    m :
        m
    n :
        n
    """
    ls = int(((n - m) / 2.0))
    cont = 0
    TCV = []
    pot = []
    a = []
    for s in range(0, (int(ls) + 1)):
        V1 = (np.power((- 1), s) * np.math.factorial((n - s)))
        V2 = ((np.math.factorial(s) * np.math.factorial((((n + m) / 2.0) - s))) * np.math.factorial((((n - m) / 2.0) - s)))
        TC = (V1 / V2)
        potencia = (n - (2.0 * s))
        TCV.append(TC)
        pot.append(potencia)
        cont = (cont + 1)
    a.append(TCV)
    a.append(pot)
    return a

def zernike_expand(L):
    """zernike_expand.

    Parameters
    ----------
    L :
        L
    """
    cont = 0
    Z = np.array([0, 0, 0, 0, 0])
    n = L
    m = L
    for i in range(0, n):
        if (cont >= L):
            break
        for j in range(0, m):
            v = ((j - i) / 2.0)
            if ((v - int(v)) == 0):
                if ((i > j) or (i == j)):
                    if (j != 0):
                        v = z_parity(cont)
                        Z = np.vstack((Z, [cont, i, j, v, np.sqrt((2.0 * (i + 1.0)))]))
                        if (cont >= L):
                            break
                        cont = (cont + 1)
                        v = z_parity(cont)
                        Z = np.vstack((Z, [cont, i, j, v, np.sqrt((2.0 * (i + 1.0)))]))
                        if (cont >= L):
                            break
                        cont = (cont + 1)
                    if (j == 0):
                        Z = np.vstack((Z, [cont, i, j, 1, np.sqrt((i + 1.0))]))
                        if (cont >= L):
                            break
                        cont = (cont + 1)
    Z = np.delete(Z, 0, 0)
    E = []
    for i in range(0, Z.shape[0]):
        (j, n, m, paR_c, raiz) = Z[i]
        a = r_zern(m, n)
        E.append(a)
    Z = Z[:L]
    return (Z, E)

def zernike_math_notation(term, Zern_pol, z_pow):
    """zernike_math_notation.

    Parameters
    ----------
    term :
        term
    Zern_pol :
        Zern_pol
    z_pow :
        z_pow
    """
    (j, n, m, par, raiz) = Zern_pol[term]
    ZZZ = ['Piston', 'Tilt x, (about y axis)', 'Tilt y, (about x axis)', 'Power or Focus', 'Astigmatism y, (45deg)', 'Astigmatism x, (0deg)', 'Coma y', 'Coma x', 'Trefoil y', 'Trefoil x', 'Primary Spherical', 'Secondary Astigmatism x', 'Secondary Astigmatism y', 'Tetrafoil x', 'Tetrafoil y', 'Secondary Coma x', 'Secondary Coma y', 'Secondary Trefoil x', 'Secondary Trefoil y', 'Pentafoil x', 'Pentafoil y', 'Secondary Spherical', 'Tertiary Astigmatism y', 'Tertiary Astigmatism x', 'Secondary Tetrafoil y', 'Secondary Tetrafoil x', '14^(1/2) (p^6) * SIN (6A)', '14^(1/2) (p^6) * COS (6A)', 'Tertiary Coma y', 'Tertiary Coma x', 'Tertiary Trefoil y', 'Tertiary Trefoil x']
    ct = z_pow[term][0]
    pot = z_pow[term][1]
    NR = ' '
    L = len(ct)
    for i in range(0, L):
        if (pot[i] == 0):
            NR = ((str(ct[i]) + '+') + NR)
        else:
            NR = ((((str(ct[i]) + 'r^') + str(int(pot[i]))) + '+') + NR)
    i = len(NR)
    NR = NR[:(- 2)]
    if (m == 1):
        mm = ''
    else:
        mm = str(int(m))
    if (par == 1):
        S = (((str(int((0.01 + (raiz * raiz)))) + '^(1/2)(') + NR) + ')')
    if (par == 3):
        S = ((((((str(int((0.01 + (raiz * raiz)))) + '^(1/2)(') + NR) + ')') + 'cos(') + mm) + 'T)')
    if (par == 2):
        S = ((((((str(int((0.01 + (raiz * raiz)))) + '^(1/2)(') + NR) + ')') + 'sin(') + mm) + 'T)')
    if (term < len(ZZZ)):
        x = ZZZ[[term][0]][:]
    else:
        x = ''
    return ((S + '   ') + x)

def Wavefront_Zernike_Phase(x, y, COEF):
    """Wavefront_Zernike_Phase.

    Parameters
    ----------
    x :
        x
    y :
        y
    COEF :
        COEF
    """
    NC = len(COEF)
    (Zern_pol, z_pow) = zernike_expand(NC)
    tcoef = COEF.shape[0]
    p = np.sqrt(((x * x) + (y * y)))
    f = np.arctan2(x, y)
    ZFP = 0.0
    for i in range(0, tcoef):
        if (COEF[i] != 0):
            ZFP = (ZFP + (COEF[i] * zernike_polynomials(i, p, f, Zern_pol, z_pow)))
    return ZFP
