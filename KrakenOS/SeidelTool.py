
import numpy as np
import KrakenOS as Kos
from .PupilTool import PupilCalc, SolveVectCross

class Seidel():
    """Seidel.
    """
    def __init__(self, PUPIL):
        """__init__.

        Parameters
        ----------
        PUPIL :
            PUPIL
        """
        self.PUPIL = PUPIL
        self.SYSTEM = self.PUPIL.SYSTEM
        self.calculate()

    def calculate(self):
        """calculate.
        """
        sup = self.PUPIL.Surf
        W = self.PUPIL.W
        ApType = self.PUPIL.ApertureType
        ApVal = self.PUPIL.ApertureValue
        fx = self.PUPIL.FieldX
        fy = self.PUPIL.FieldY
        field = np.sqrt(((fx ** 2) + (fy ** 2)))
        fieldType = self.PUPIL.FieldType
        Pup = PupilCalc(self.SYSTEM, sup, W, ApType, ApVal)
        Prx = self.SYSTEM.Parax(W)
        (SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP, C, NN, d) = Prx
        n_1 = np.asarray(NN)
        n_sign = np.sign(NN)
        n_abs = np.abs(NN)
        for si in range(0, len(n_sign)):
            n_abs[si:] = (n_sign[si] * n_abs[si:])
        n_1 = np.copy(n_abs)
        n_2 = np.copy(n_abs)
        n_1 = np.insert(n_1, 0, n_1[0])
        n_elements = self.SYSTEM.n
        c = []
        d = []
        kon = []
        asp2 = []
        asp4 = []
        asp6 = []
        d.append(0)
        for i in range(0, n_elements):
            rc = self.SYSTEM.SDT[i].Rc
            if (rc == 0):
                c.append(rc)
            else:
                c.append((1 / rc))
            t = self.SYSTEM.SDT[i].Thickness
            d.append(t)
            kon.append(self.SYSTEM.SDT[i].k)
            asp2.append(self.SYSTEM.SDT[i].AspherData[0])
            asp4.append(self.SYSTEM.SDT[i].AspherData[1])
            asp6.append(self.SYSTEM.SDT[i].AspherData[2])
        c.append(0)
        d.append(0)
        kon.append(0)
        asp2.append(0)
        asp4.append(0)
        asp6.append(0)
        t = d[1]
        d[1] = (d[1] - Pup.PosPupInp[2])
        Diam = (Pup.RadPupInp * 2.0)
        c = np.asarray(c)
        epsilon = np.asarray(kon)
        asp2 = np.asarray(asp2)
        asp4 = np.asarray(asp4)
        asp6 = np.asarray(asp6)
        c_orig = np.copy(c)
        c = (c + (2.0 * asp2))
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
        if (fieldType == 'angle'):
            Theta_marg = 0
            Theta_bar = (- np.deg2rad(field))
        PosObj = (t + d[1])
        if (fieldType == 'height'):
            Theta_marg = np.arctan2(Pup.RadPupInp, PosObj)
            Theta_bar = (- np.arctan2(field, PosObj))
        u.append(Theta_marg)
        h.append((Diam / 2.0))
        u_bar.append(Theta_bar)
        h_bar.append(0)
        E.append(0)
        for k in range(0, (len(c) - 1)):
            nk = n_1[k]
            nkm1 = n_2[k]
            u.append((((nk * u[k]) - ((h[k] * c[k]) * (nkm1 - nk))) / nkm1))
            h.append((h[k] + (u[(k + 1)] * d[(k + 1)])))
            u_bar.append((((nk * u_bar[k]) - ((h_bar[k] * c[k]) * (nkm1 - nk))) / nkm1))
            h_bar.append((h_bar[k] + (u_bar[(k + 1)] * d[(k + 1)])))
            Epro = (((- d[(k + 1)]) / ((nkm1 * h[(k + 1)]) * h[k])) + E[k])
            E.append(Epro)
            Delta_unn.append(((u[(k + 1)] / nkm1) - (u[k] / nk)))
            Delta_nn.append(((1.0 / nkm1) - (1.0 / nk)))
            Delta_n.append((nkm1 - nk))
            Delta_1sn.append(((1.0 / (nkm1 ** 2)) - (1 / (nk ** 2))))
            A.append((nk * ((h[k] * c[k]) + u[k])))
        Abar = (1.0 * ((h_bar[0] * c[0]) + u_bar[0]))
        H = ((A[0] * h_bar[0]) - (Abar * h[0]))
        for k in range(0, (len(c) - 1)):
            A_bar.append(((H / h[k]) * (((A[k] * h[k]) * E[k]) - 1)))
            sI.append((((A[k] ** 2.0) * h[k]) * Delta_unn[k]))
            sII.append((((A[k] * A_bar[k]) * h[k]) * Delta_unn[k]))
            sIII.append((((A_bar[k] ** 2.0) * h[k]) * Delta_unn[k]))
            sIV.append((((H ** 2.0) * c[k]) * Delta_nn[k]))
            if (A[k] != 0):
                P = (c[k] * Delta_nn[k])
                p1 = ((H ** 2.0) * P)
                p2 = (((A_bar[k] ** 2.0) * h[k]) * Delta_unn[k])
                sV.append(((A_bar[k] / A[k]) * (p1 + p2)))
            else:
                sV.append(0)
            a = ((epsilon[k] * (c[k] ** 3.0)) * (h[k] ** 4.0))
            sI_k.append(((a * Delta_n[k]) * ((h_bar[k] / h[k]) ** 0)))
            sII_k.append(((a * Delta_n[k]) * ((h_bar[k] / h[k]) ** 1)))
            sIII_k.append(((a * Delta_n[k]) * ((h_bar[k] / h[k]) ** 2)))
            sIV_k.append((a * 0.0))
            sV_k.append(((a * Delta_n[k]) * ((h_bar[k] / h[k]) ** 3)))
        sI_a = []
        sII_a = []
        sIII_a = []
        sIV_a = []
        sV_a = []
        c = c_orig
        Asp = (asp4 - ((asp2 / 4.0) * (((4.0 * (asp2 ** 2)) + ((6.0 * asp2) * c)) + (3.0 * (c ** 2.0)))))
        for k in range(0, (len(c) - 1)):
            a = (((8.0 * Asp[k]) * (h[k] ** 4.0)) * Delta_n[k])
            sI_a.append((a * ((h_bar[k] / h[k]) ** 0)))
            sII_a.append((a * ((h_bar[k] / h[k]) ** 1)))
            sIII_a.append((a * ((h_bar[k] / h[k]) ** 2)))
            sIV_a.append((a * 0.0))
            sV_a.append((a * ((h_bar[k] / h[k]) ** 3)))
        sI = np.asarray(sI)
        sII = np.asarray(sII)
        sIII = np.asarray(sIII)
        sIV = np.asarray(sIV)
        sV = np.asarray(sV)
        sI_k = np.asarray(sI_k)
        sII_k = np.asarray(sII_k)
        sIII_k = np.asarray(sIII_k)
        sIV_k = np.asarray(sIV_k)
        sV_k = np.asarray(sV_k)
        sI_a = np.asarray(sI_a)
        sII_a = np.asarray(sII_a)
        sIII_a = np.asarray(sIII_a)

        self.si = ((sI - sI_k) - sI_a)
        self.sii = ((sII - sII_k) - sII_a)
        self.siii = ((sIII - sIII_k) - sIII_a)
        self.siv = ((sIV - sIV_k) - sIV_a)
        self.sv = ((sV - sV_k) - sV_a)
        #---------------------------------#
        self.SAC = [self.si, self.sii, self.siii, self.siv, self.sv]
        self.SAC_TOTAL = [np.sum(self.si), np.sum(self.sii), np.sum(self.siii), np.sum(self.siv), np.sum(self.sv)]
        self.SAC_NM = ['si', 'sii', 'siii', 'siv', 'sv']

        Units = 1000.0
        self.W040 = (((Units * (1 / 8.0)) * self.si) / W)
        self.W131 = (((Units * (1 / 2.0)) * self.sii) / W)
        self.W222 = (((Units * (1 / 2.0)) * self.siii) / W)
        self.W220 = (((Units * (1 / 4.0)) * self.siv) / W)
        self.W311 = (((Units * (1 / 2.0)) * self.sv) / W)
        #---------------------------------#
        self.SCW = [self.W040, self.W131, self.W222, self.W220, self.W311]
        self.SCW_TOTAL = [np.sum(self.W040), np.sum(self.W131), np.sum(self.W222), np.sum(self.W220), np.sum(self.W311)]
        self.SCW_NM = ['W040', 'W131', 'W222', 'W220', 'W311']

        self.TSPH = (self.si / ((2.0 * u[(- 1)]) * n_1[(- 1)]))
        self.TSCO = (self.sii / ((2.0 * u[(- 1)]) * n_1[(- 1)]))
        self.TTCO = ((3.0 * self.sii) / ((2.0 * u[(- 1)]) * n_1[(- 1)]))
        self.TAST = ((- self.siii) / (u[(- 1)] * n_1[(- 1)]))
        self.TPFC = ((- self.siv) / ((2.0 * u[(- 1)]) * n_1[(- 1)]))
        self.TSFC = ((self.siii + self.siv) / ((2.0 * u[(- 1)]) * n_1[(- 1)]))
        self.TTFC = (((3.0 * self.siii) + self.siv) / ((2.0 * u[(- 1)]) * n_1[(- 1)]))
        self.TDIS = ((- self.sv) / ((2.0 * u[(- 1)]) * n_1[(- 1)]))
        #---------------------------------#

        self.TAC = [self.TSPH, self.TSCO, self.TTCO, self.TAST, self.TPFC, self.TSFC, self.TTFC, self.TDIS]
        self.TAC_TOTAL = [np.sum(self.TSPH), np.sum(self.TSCO), np.sum(self.TTCO), np.sum(self.TAST), np.sum(self.TPFC), np.sum(self.TSFC), np.sum(self.TTFC), np.sum(self.TDIS)]
        self.TAC_NM = ['TSPH', 'TSCO', 'TTCO', 'TAST', 'TPFC', 'TSFC', 'TTFC', 'TDIS']

        #---------------------------------#
        self.LSPH = (self.si / (((2.0 * u[(- 1)]) * u[(- 1)]) * n_1[(- 1)]))
        self.LSCO = (self.sii / (((2.0 * u[(- 1)]) * u[(- 1)]) * n_1[(- 1)]))
        self.LTCO = ((3.0 * self.sii) / (((2.0 * u[(- 1)]) * u[(- 1)]) * n_1[(- 1)]))
        self.LAST = (self.siii / ((u[(- 1)] * u[(- 1)]) * n_1[(- 1)]))
        self.LPFC = (self.siv / (((2.0 * u[(- 1)]) * u[(- 1)]) * n_1[(- 1)]))
        self.LSFC = ((self.siii + self.siv) / (((2.0 * u[(- 1)]) * u[(- 1)]) * n_1[(- 1)]))
        self.LTFC = (((3.0 * self.siii) + self.siv) / (((2.0 * u[(- 1)]) * u[(- 1)]) * n_1[(- 1)]))
        self.LDIS = ((- self.sv) / (((2.0 * u[(- 1)]) * u[(- 1)]) * n_1[(- 1)]))
        #---------------------------------#
        self.LAC = [self.LSPH, self.LSCO, self.LTCO, self.LAST, self.LPFC, self.LSFC, self.LTFC, self.LDIS]
        self.LAC_TOTAL = [np.sum(self.LSPH), np.sum(self.LSCO), np.sum(self.LTCO), np.sum(self.LAST), np.sum(self.LPFC), np.sum(self.LSFC), np.sum(self.LTFC), np.sum(self.LDIS)]
        self.LAC_NM = ['LSPH', 'LSCO', 'LTCO', 'LAST', 'LPFC', 'LSFC', 'LTFC', 'LDIS']

        self.SAC_AN = 'Seidel Aberration Coefficents'
        self.SCW_AN = 'Seidel coefficients in waves'
        self.TAC_AN = 'Transverse Aberration Coefficents'
        self.LAC_AN = 'Longitudinal Aberration Coefficents'
        self.AberNames = [self.SAC_AN, self.SCW_AN, self.TAC_AN, self.LAC_AN]

