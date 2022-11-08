
import numpy as np
import math
from .SurfTools import surface_tools as SUT

class InterNormalCalc():
    """InterNormalCalc.
    """

    def __init__(self, SurfData, TypTot, PR3D, HS):
        """__init__.

        Parameters
        ----------
        SurfData :
            SurfData
        TypTot :
            TypTot
        PR3D :
            PR3D
        HS :
            HS
        """
        self.HS = HS
        self.SDT = SurfData
        self.n = len(self.SDT)
        self.SuTo = SUT(SurfData)
        self.Pr3D = PR3D
        self.Disable_Inner = 1
        self.ExtraDiameter = 0
        self.AAA = self.Pr3D.AAA
        self.BBB = self.Pr3D.BBB
        self.DDD = self.Pr3D.DDD
        self.EEE = self.Pr3D.EEE
        self.GlassOnSide = self.Pr3D.GlassOnSide
        self.side_number = self.Pr3D.side_number
        self.TypeTotal = self.Pr3D.TypeTotal
        self.TRANS_1A = self.Pr3D.TRANS_1A
        self.TRANS_2A = self.Pr3D.TRANS_2A
        self.Pn = np.asarray([0.0, 0.0, 0.0])

        self.P1 = np.asarray([0.0, 0.0, 0.0, 1.0])
        self.P2 = np.asarray([0.0, 0.0, 0.0, 1.0])
        self.P_z1 = 10000000.0

    def __SigmaHitTransfSpace(self, PP_start, PP_stop, j):
        """__SigmaHitTransfSpace.

        Parameters
        ----------
        PP_start :
            PP_start
        PP_stop :
            PP_stop
        j :
            j
        """

        StopPoint = np.array([PP_stop[0], PP_stop[1], PP_stop[2], 1.0])
        StarPoint = np.array([PP_start[0], PP_start[1], PP_start[2], 1.0])

        SurfHit = 1
        P_SurfHit = self.Pr3D.TRANS_1A[j].dot(StopPoint)
        Px1 = P_SurfHit[(0, 0)]
        Py1 = P_SurfHit[(0, 1)]
        Pz1 = P_SurfHit[(0, 2)]

        P_start = self.Pr3D.TRANS_1A[j].dot(StarPoint)
        P_x1 = P_start[(0, 0)]
        P_y1 = P_start[(0, 1)]
        P_z1 = P_start[(0, 2)]

        P12 = [(Px1 - P_x1), (Py1 - P_y1), (Pz1 - P_z1)]
        [L, M, N] = (P12 / np.linalg.norm(P12))

        Px1 = (((L / N) * (- P_z1)) + P_x1)
        Py1 = (((M / N) * (- P_z1)) + P_y1)

        Pz1 = 0
        SurfHit = 1

        P_x2 = 0
        P_y2 = 0
        P_z2 = 0

        if (len(np.shape(P_x1)) != 0):
            L = L[0]
            M = M[0]
            N = N[0]
            P_x1 = P_x1[0]
            P_y1 = P_y1[0]
            P_z1 = P_z1[0]
            Px1 = Px1[0]
            Py1 = Py1[0]
            Pz1 = Pz1[0]
            P_x2 = P_x2[0]
            P_y2 = P_y2[0]
            P_z2 = P_z2[0]
            SurfHit = SurfHit[0]

        if (self.SDT[j].Thin_Lens == 0):
            self.vj = j

            ASD=np.sqrt(((Px1-self.SDT[j].SubAperture[2])**2) + ((Py1-self.SDT[j].SubAperture[1])**2))
            D0 = (2.0 * ASD)
            DiamInf = ((self.SDT[j].InDiameter * self.SDT[j].SubAperture[0]) * self.Disable_Inner)
            DiamSup = ((self.SDT[j].Diameter * self.SDT[j].SubAperture[0]) + (10000.0 * self.ExtraDiameter))


            if ((D0 > DiamSup) or (D0 < DiamInf)):
                SurfHit = 0
                P_x2 = 0
                P_y2 = 0
                P_z2 = 0
            else:
                (P_x2, P_y2, P_z2) = self.HS.SolveHit(Px1, Py1, Pz1, L, M, N, j)
                if (not math.isnan(P_z2)):
                    P_x2 = self.HS.vevaX
                    P_y2 = self.HS.vevaY

                    ASD=np.sqrt(((Px1-self.SDT[j].SubAperture[2])**2) + ((Py1-self.SDT[j].SubAperture[1])**2))
                    D0 = (2.0 * ASD)

                    DiamInf = ((self.SDT[j].InDiameter * self.SDT[j].SubAperture[0]) * self.Disable_Inner)
                    DiamSup = ((self.SDT[j].Diameter * self.SDT[j].SubAperture[0]) + (10000.0 * self.ExtraDiameter))
                    if ((D0 > DiamSup) or (D0 < DiamInf)):
                        SurfHit = 0
                else:
                    SurfHit = 0
                    P_x2 = 0
                    P_y2 = 0
                    P_z2 = 0
        else:
            ASD=np.sqrt(((Px1 - self.SDT[j].SubAperture[2])**2) + ((Py1 - self.SDT[j].SubAperture[1])**2))
            D0 = (2.0 * ASD)*0.999999

            if ((D0 > self.SDT[j].Diameter * self.SDT[j].SubAperture[0]) or (D0 < self.SDT[j].InDiameter * self.SDT[j].SubAperture[0] )):
                SurfHit = 0
                P_x2 = 0
                P_y2 = 0
                P_z2 = 0

            else:
                P_x2 = ((L / N) * self.SDT[j].Thin_Lens)
                P_y2 = ((M / N) * self.SDT[j].Thin_Lens)
                P_z2 = self.SDT[j].Thin_Lens

        return (SurfHit, P_x2, P_y2, P_z2, Px1, Py1, Pz1, L, M, N)

    def __SigmaHitTransfSpaceFast(self, PP_start, PP_stop, j):

        StopPoint = np.array([PP_stop[0], PP_stop[1], PP_stop[2], 1.0])
        StarPoint = np.array([PP_start[0], PP_start[1], PP_start[2], 1.0])

        SurfHit = 1
        P_SurfHit = self.Pr3D.TRANS_1A[j].dot(StopPoint)
        Px1 = P_SurfHit[(0, 0)]
        Py1 = P_SurfHit[(0, 1)]
        Pz1 = P_SurfHit[(0, 2)]

        P_start = self.Pr3D.TRANS_1A[j].dot(StarPoint)
        P_x1 = P_start[(0, 0)]
        P_y1 = P_start[(0, 1)]
        P_z1 = P_start[(0, 2)]

        P12 = [(Px1 - P_x1), (Py1 - P_y1), (Pz1 - P_z1)]
        [L, M, N] = (P12 / np.linalg.norm(P12))

        Px1 = (((L / N) * (- P_z1)) + P_x1)
        Py1 = (((M / N) * (- P_z1)) + P_y1)

        Pz1 = 0
        SurfHit = 1

        P_x2 = 0
        P_y2 = 0
        P_z2 = 0


        (P_x2, P_y2, P_z2) = self.HS.SolveHit(Px1, Py1, Pz1, L, M, N, j)
        if (not math.isnan(P_z2)):
            P_x2 = self.HS.vevaX
            P_y2 = self.HS.vevaY

        else:
            SurfHit = 0
            P_x2 = 0
            P_y2 = 0
            P_z2 = 0

        return (SurfHit, P_x2, P_y2, P_z2, Px1, Py1, Pz1)

    def __ParaxCalcObjOut2OrigSpace(self, Px2, Py2, Pz2, Px1, Py1, Pz1, j):
        """__ParaxCalcObjOut2OrigSpace.

        Parameters
        ----------
        Px2 :
            Px2
        Py2 :
            Py2
        Pz2 :
            Pz2
        Px1 :
            Px1
        Py1 :
            Py1
        Pz1 :
            Pz1
        j :
            j
        """
        P1 = [Px1, Py1, Pz1, 1]
        P2 = [Px2, Py2, Pz2, 1]
        NP1 = self.TRANS_2A[j].dot(P1)
        NP2 = self.TRANS_2A[j].dot(P2)
        Pn = np.asarray([(- (NP1[(0, 0)] - NP2[(0, 0)])), (- (NP1[(0, 1)] - NP2[(0, 1)])), (- (NP1[(0, 2)] - NP2[(0, 2)]))])
        norm = (Pn / np.linalg.norm(Pn))
        PTO_exit = [NP1[(0, 0)], NP1[(0, 1)], NP1[(0, 2)]]
        PTO_exit_Object_Space = [Px1, Py1, Pz1]
        return (norm, PTO_exit, PTO_exit_Object_Space)

    def __SigmaOutOrigSpace(self, P_x2, P_y2, P_z2, j):
        """__SigmaOutOrigSpace.

        Parameters
        ----------
        P_x2 :
            P_x2
        P_y2 :
            P_y2
        P_z2 :
            P_z2
        j :
            j
        """
        (New_L, New_M, New_N) = self.HS.SurfDer(P_x2, P_y2, P_z2)

        Pz1z2 = (self.P_z1 - P_z2)

        P_x1 = ((Pz1z2 * (New_L / New_N)) + P_x2)
        P_y1 = ((Pz1z2 * (New_M / New_N)) + P_y2)


        self.P1[0], self.P1[1], self.P1[2] = P_x1, P_y1, self.P_z1
        self.P2[0], self.P2[1], self.P2[2] = P_x2, P_y2, P_z2

        NP1 = self.TRANS_2A[j].dot(self.P1)
        NP2 = self.TRANS_2A[j].dot(self.P2)


        self.Pn[0] = - (NP1[(0, 0)] - NP2[(0, 0)])
        self.Pn[1] = - (NP1[(0, 1)] - NP2[(0, 1)])
        self.Pn[2] = - (NP1[(0, 2)] - NP2[(0, 2)])


        LNOR=np.sqrt((self.Pn[0]**2.)+(self.Pn[1]**2.)+(self.Pn[2]**2.))

        norm = (self.Pn / LNOR)

        PTO_exit = [NP2[(0, 0)], NP2[(0, 1)], NP2[(0, 2)]]
        PTO_exit_Object_Space = [P_x2, P_y2, P_z2]


        return (norm, PTO_exit, PTO_exit_Object_Space)

    def __HitOnMask(self, PP_start, PP_stop, j):
        """__HitOnMask.

        Parameters
        ----------
        PP_start :
            PP_start
        PP_stop :
            PP_stop
        j :
            j
        """
        SurfHit = 1
        HITS_CONT = []
        if (self.SDT[j].Mask_Type != 0):
            OBJECT = self.DDD[j]
            for obj in OBJECT:
                (inter_mask, ind_mask) = obj.ray_trace(PP_start, PP_stop)
                Hit_MASK = np.shape(inter_mask)[0]
                HITS_CONT.append(Hit_MASK)
            HITS_CONT = np.asarray(HITS_CONT)
            if np.any((HITS_CONT == 1)):
                SurfHit_MASK = 1
            else:
                SurfHit_MASK = 0
            if (self.SDT[j].Mask_Type == 1):
                if (SurfHit_MASK == 1):
                    SurfHit = 1
                else:
                    SurfHit = 0
            if (self.SDT[j].Mask_Type == 2):
                if (SurfHit_MASK == 1):
                    SurfHit = 0
                else:
                    SurfHit = 1
        return SurfHit

    def __GrooveDirectionVector(self, j):
        """__GrooveDirectionVector.

        Parameters
        ----------
        j :
            j
        """
        self.P1[0], self.P1[1], self.P1[2] = 0, 0, 0
        self.P2[0], self.P2[1], self.P2[2] = -np.cos(np.deg2rad(self.SDT[j].Grating_Angle)), -np.sin(np.deg2rad(self.SDT[j].Grating_Angle)),0

        NP1 = self.TRANS_2A[j].dot(self.P1)
        NP2 = self.TRANS_2A[j].dot(self.P2)

        self.Pn[0] = - (NP1[(0, 0)] - NP2[(0, 0)])
        self.Pn[1] = - (NP1[(0, 1)] - NP2[(0, 1)])
        self.Pn[2] = - (NP1[(0, 2)] - NP2[(0, 2)])

        LNOR=np.sqrt((self.Pn[0]**2.)+(self.Pn[1]**2.)+(self.Pn[2]**2.))

        Pg_v = (self.Pn / LNOR)
        return Pg_v

    def InterNormal(self, PP_start, PP_stop, j, jj):
        """InterNormal.

        Parameters
        ----------
        PP_start :
            PP_start
        PP_stop :
            PP_stop
        j :
            j
        jj :
            jj
        """
        PTO_exit = [0, 0, 0]
        PTO_exit_Object_Space = [0, 0, 0]
        LMN_exit_Object_Space = [0, 0, 1]
        norm = [0, 0, 1]
        SurfHit = 1

        if (self.SDT[j].Diff_Ord == 0):
            Pgn = [0, 1, 0]

        else:
            Pgn = self.__GrooveDirectionVector(j)

        if (self.TypeTotal[jj] == 0):
            SurfHit = 1
            SurfHit = self.__HitOnMask(PP_start, PP_stop, j)

            if (SurfHit != 0):

                (SurfHit, Px2, Py2, Pz2, Px1, Py1, Pz1, L, M, N) = self.__SigmaHitTransfSpace(PP_start, PP_stop, j)
                LMN_exit_Object_Space = [L, M, N]

                if (self.SDT[j].Thin_Lens == 0):
                    (norm, PTO_exit, PTO_exit_Object_Space) = self.__SigmaOutOrigSpace(Px2, Py2, Pz2, j)

                else:
                    (norm, PTO_exit, PTO_exit_Object_Space) = self.__ParaxCalcObjOut2OrigSpace(Px2, Py2, Pz2, Px1, Py1, Pz1, j)
        else:
            (SurfHit, norm, PTO_exit, Pgn) = self.__InterNormalSolidObject(jj, PP_start, PP_stop)
        return (SurfHit, np.asarray(norm), np.asarray(PTO_exit), np.asarray(Pgn), np.asarray(PTO_exit_Object_Space), np.asarray(LMN_exit_Object_Space), j)

    def __InterNormalSolidObject(self, jj, PP_start, PP_stop):
        """__InterNormalSolidObject.

        Parameters
        ----------
        jj :
            jj
        PP_start :
            PP_start
        PP_stop :
            PP_stop
        """
        Pgn = np.asarray([0, 1, 0])
        PTO_exit = [0, 0, 0]
        norm = [0, 0, 1]
        (inter, ind) = self.EEE[jj].ray_trace(PP_start, PP_stop)
        SurfHit = np.shape(inter)[0]
        if (SurfHit != 0):
            s = 0
            h = []
            for f in ind:
                PD = (np.asarray(inter[s]) - np.asarray(PP_start))
                distance = np.linalg.norm(PD)
                if (np.abs(distance) < 0.05):
                    distance = 99999999999999.9
                h.append(distance)
                s = (s + 1)
            index = np.argmin(np.asarray(h))
            PTO_exit = inter[index]
            NOR = self.EEE[jj].cell_normals
            norm = NOR[ind[index]]
            Pgn = np.asarray([0, 0, 1])
        return (SurfHit, norm, PTO_exit, Pgn)
