# from numba import jit

import numpy as np
import os
import sys
import random
import inspect

from .SurfTools import surface_tools as SUT
from .Prerequisites3D import *

from .Physics import *
from .HitOnSurf import *
from .InterNormalCalc import *

import timeit



# from .__init__ import RUTE as ModRute

def prob(pro):
    """prob.

    Parameters
    ----------
    pro :
        pro
    """
    a_list = [0, 1]
    prob = pro
    distribution = [prob, (1.0 - prob)]
    random_number = random.choices(a_list, distribution)
    return random_number

class system():


    # print(ModRute)

    """system.
            SYSTEM CLASS ATRIBUTES AND IMPLEMENTATIONS:


          system.Trace (pS, dC, wV)
             Sequential ray tracing.
             pS = [1.0, 0.0, 0.0] – Ray origin coordinates
             dC = [0.0,0.0,1.0] - The directing cosines
             wV = 0.4 - Wavelength


          system.NsTrace (pS, dC, wV)
             Non-Sequential ray tracing


          system.Parax (w)
             Paraxial optics calculations


          system.disable_inner
             Enables the central aperture.


          system.enable_inner
             Disables the central aperture.


          system.SURFACE
	    Returns the surfaces the ray passed through.


          system.NAME
        	    Returns surface names that the ray passed through


          system.GLASS
        	    Returns materials that the ray passed through.


          system.XYZ
             [X, Y, Z] ray coordinates from its origin to the image plane.


          system.OST_XYZ
             [X, Y, Z] coordinates of ray intersections in reference to
            a coordinate system at its vertex, even if this vertex has
            a translation or rotation.


          system.DISTANCE
             List of distances traveled by the ray.


          system.OP
        	List of optical paths.


          system.TOP
        	Total optical path.


          system.TOP_S     List of the ray's optical path by sections.


          system.ALPHA
             List the materials absorption coefficients


          system.BULK_TRANS
             List the transmission through all the system absorption
             coefficients are considered.


          system.S_LMN
             Surface normal direction cosines [L, M, N].


          system.LMN
             incident ray direction cosines [L, M, N].


          system.R_LMN
             Resulting ray direction cosines [L, M, N].


          system.N0
             Refractive indices before and after each interface


          system.N1
             Refractive indices after each interface.
             This is useful to differentiate between index
             before and after an iteration. Example:

             N0 = [n1, n2, n3, n4, n5]
             N1 = [n2, n3, n4, n5, n5]


          system.WAV
        	    Wavelength of the ray (µm)


          system.G_LMN
             [L, M, N] Direction cosines that define the lines
             on the diffraction grating on the plane.


          system.ORDER
        	     Ray diffraction order.


          system.GRATING_D
	     Distance between lines of the diffraction grating.
             Units (Microns)


          system.RP
	     Fresnel reflection coefficients for polarization P.


          system.RS
             Fresnel reflection coefficients for polarization S.


          system.TP
            Fresnel transmission coefficients for polarization P.


          system.TS
            Fresnel transmission coefficients for polarization S.


          system.TTBE
            Total energy transmitted or reflected per element.


          system.TT
            Total energy transmitted or reflected total.


          system.targ_surf (int)
            Limits the ray tracing to the defined surface


          system.flat_surf (int)
            Change a surface to flat.

    """


    def __init__(self, SurfData, KN_Setup):
        """__init__.

        Parameters
        ----------
        SurfData :
            SurfData
        KN_Setup :
            KN_Setup
        """
        self.ExectTime=[]

        self.SDT = SurfData
        self.update = False
        self.S_Matrix = []
        self.N_Matrix = []
        self.SistemMatrix = np.matrix([[1.0, 0.0], [0.0, 1.0]])
        (self.a, self.b, self.c, self.d, self.EFFL, self.PPA, self.PPP) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.SETUP = KN_Setup
        self.n = len(self.SDT)
        for ii in range(0, self.n):
            self.SDT[ii].SaveSetup()
        self.SuTo = SUT(self.SDT)
        self.Object_Num = np.arange(0, self.n, 1)
        self.Targ_Surf = self.n
        self.SuTo.Surface_Flattener = 0
        self.Disable_Inner = 1
        self.ExtraDiameter = 0
        self.SuTo.ErrSurfCase = 0
        self.__SurFuncSuscrip()
        self.Pr3D = Prerequisites(self.SDT, self.SuTo)
        self.__PrerequisitesGlass()
        self.Pr3D.Prerequisites3SMath()
        self.Pr3D.Prerequisites3DSolids()
        self.PreWave = (- 1e-09)
        self.AAA = self.Pr3D.AAA
        self.BBB = self.Pr3D.BBB
        self.DDD = self.Pr3D.DDD
        self.EEE = self.Pr3D.EEE
        self.GlassOnSide = self.Pr3D.GlassOnSide
        self.side_number = self.Pr3D.side_number
        self.TypeTotal = self.Pr3D.TypeTotal
        self.TRANS_1A = self.Pr3D.TRANS_1A
        self.TRANS_2A = self.Pr3D.TRANS_2A
        self.HS = Hit_Solver(self.SDT)
        self.INORM = InterNormalCalc(self.SDT, self.TypeTotal, self.Pr3D, self.HS)
        self.INORM.Disable_Inner = 1
        self.Pr3D.Disable_Inner = 1
        (self.c_p, self.n_p, self.d_p) = (0, 0, 0)
        self.tt=1.

        # (ResVec, CurrN, sign) = self.SDT[2].PHYSICS.calculate([-0.01295049, 0.14862896 ,0.98880823], [-0., -0., -1.], 1.0, 1.0, [0., 1. ,0.], 0.0, 0.0, 0.6, 0.)

    def __SurFuncSuscrip(self):
        """__SurFuncSuscrip.
        """
        for i in range(0, self.n):
            self.SDT[i].build_surface_function()

    def __PrerequisitesGlass(self):
        """__PrerequisitesGlass.
        """
        self.Glass = []
        self.GlobGlass = []
        for i in range(0, self.n):
            if (type(self.SDT[i].Glass) == type(1.0)):
                self.SDT[i].Glass = ('AIR_' + str(self.SDT[i].Glass))
                print(self.SDT[i].Glass)
            self.Glass.append(self.SDT[i].Glass.replace(' ', ''))
            self.GlobGlass.append(self.SDT[i].Glass.replace(' ', ''))
            if (self.GlobGlass[i] == 'NULL'):
                self.GlobGlass[i] = self.GlobGlass[(i - 1)]

    def __NonSequentialChooserToot(self, A_RayOrig, A_Proto_pTarget, k):
        """__NonSequentialChooserToot.

        Parameters
        ----------
        A_RayOrig :
            A_RayOrig
        A_Proto_pTarget :
            A_Proto_pTarget
        k :
            k
        """
        ng = self.GlassOnSide[k]
        A_Glass = self.SDT[ng].Glass
        if (A_Glass == 'NULL'):
            distance = 99999999999999.9
        else:
            (A_pTarget, A_SurfHit) = self.EEE[k].ray_trace(A_RayOrig, A_Proto_pTarget)
            if (len(A_SurfHit) == 0):
                distance = 99999999999999.9
            else:
                s = 0
                h = []
                for f in A_SurfHit:
                    PD = (np.asarray(A_pTarget[s]) - np.asarray(A_RayOrig))
                    distance = np.linalg.norm(PD)
                    if (np.abs(distance) < 0.05):
                        distance = 99999999999999.9
                    h.append(distance)
                    s = (s + 1)
                distance = np.min(np.asarray(h))
        return distance

    def __NonSequentialChooser(self, SIGN, A_RayOrig, ResVec, j):
        """__NonSequentialChooser.

        Parameters
        ----------
        SIGN :
            SIGN
        A_RayOrig :
            A_RayOrig
        ResVec :
            ResVec
        j :
            j
        """
        chooser = []
        [SLL, SMM, SNN] = ResVec
        distance = 99999999999999.9
        [Px, Py, Pz] = A_RayOrig
        [LL, MM, NN] = ResVec
        A_Proto_pTarget = (np.asarray(A_RayOrig) + ((np.asarray(ResVec) * 999999999.9) * SIGN))
        for k in range(1, len(self.EEE)):
            distance = self.__NonSequentialChooserToot(A_RayOrig, A_Proto_pTarget, k)
            chooser.append(distance)
        chooser = np.asarray(chooser)
        jj = (np.argmin(chooser) + 1)
        chooser[(jj - 1)] = 99999999999999.9
        kk = (np.argmin(chooser) + 1)
        (A_pTarget, A_SurfHit) = self.EEE[int(jj)].ray_trace(A_RayOrig, A_Proto_pTarget)
        PRR = np.shape(A_SurfHit)[0]
        return (int(j), int(jj), int(kk), PRR)

    def __CollectDataInit(self):
        """__CollectDataInit.
        """
        self.val = 1
        self.SURFACE = []
        self.NAME = []
        self.GLASS = []
        self.S_XYZ = []
        self.T_XYZ = []
        self.XYZ = []
        self.XYZ.append([0, 0, 0])
        self.OST_XYZ = []
        self.S_LMN = []
        self.LMN = []
        self.R_LMN = []
        self.N0 = []
        self.N1 = []
        self.WAV = 1.0
        self.G_LMN = []
        self.ORDER = []
        self.GRATING = []
        self.DISTANCE = []
        self.OP = []
        self.TOP_S = []
        self.TOP = 0
        self.ALPHA = [0.0]
        self.BULK_TRANS = []
        self.RP = []
        self.RS = []
        self.TP = []
        self.TS = []
        self.TTBE = []
        self.TT = 1.0
        return None

    def __EmptyCollect(self, pS, dC, WaveLength, j):
        """__EmptyCollect.

        Parameters
        ----------
        pS :
            pS
        dC :
            dC
        WaveLength :
            WaveLength
        j :
            j
        """
        Empty = np.asarray([])
        RayTraceType = 0
        ValToSav = [Empty, Empty, pS, pS, Empty, Empty, dC, Empty, Empty, Empty, WaveLength, Empty, Empty, Empty, Empty, j, RayTraceType]
        self.__CollectData(ValToSav)

    def __WavePrecalc(self):
        """__WavePrecalc.
        """
        if (self.Wave != self.PreWave):
            self.N_Prec = []
            self.AlphaPrecal = []
            self.PreWave = self.Wave
            for i in range(0, self.n):
                (NP, AP) = n_wave_dispersion(self.SETUP, self.GlobGlass[i], self.Wave)
                self.N_Prec.append(NP)
                self.AlphaPrecal.append(AP)

    def __CollectData(self, ValToSav):
        """__CollectData.

        Parameters
        ----------
        ValToSav :
            ValToSav
        """
        [Glass, alpha, RayOrig, pTarget, HitObjSpace, SurfNorm, ImpVec, ResVec, PrevN, CurrN, WaveLength, D, Ord, GrSpa, Name, j, RayTraceType] = ValToSav



        self.SURFACE.append(j)
        self.NAME.append(Name)
        self.GLASS.append(Glass)
        self.S_XYZ.append(RayOrig)
        self.T_XYZ.append(pTarget)
        self.XYZ[0] = self.S_XYZ[0]
        self.XYZ.append(pTarget)
        self.OST_XYZ.append(HitObjSpace)
        p = (np.asarray(RayOrig) - np.asarray(pTarget))
        dist = np.linalg.norm(p)
        self.DISTANCE.append(dist)
        self.OP.append((dist * PrevN))

        self.TOP = (self.TOP + (dist * PrevN))
        self.TOP_S.append(self.TOP)
        self.ALPHA.append(alpha)
        self.S_LMN.append(SurfNorm)
        self.LMN.append(ImpVec)
        self.R_LMN.append(ResVec)
        self.N0.append(PrevN)
        self.N1.append(CurrN)
        self.WAV = WaveLength
        self.G_LMN.append(D)
        self.ORDER.append(Ord)
        self.GRATING.append(GrSpa)
        if (self.val == 1):
            (Rp, Rs, Tp, Ts) = FresnelEnergy(Glass, PrevN, CurrN, ImpVec, SurfNorm, ResVec, self.SETUP, self.Wave)
        else:
            (Rp, Rs, Tp, Ts) = (0, 0, 0, 0)
        self.RP.append(Rp)
        self.RS.append(Rs)
        self.TP.append(Tp)
        self.TS.append(Ts)
        if ((RayTraceType == 0) or (RayTraceType == 1)):
            if (Glass == 'MIRROR'):
                self.tt = ((1.0 * (Rp + Rs)) / 2.0)
                self.BULK_TRANS.append(self.tt)
            if (Glass != 'MIRROR'):
                IT = np.exp(((- self.ALPHA[-2]) * dist))

                self.BULK_TRANS.append(IT)
                self.tt = (Tp + Ts) / 2.0
        else:
            self.tt = 1.0


        # si self.tt está vacio entonces es cero #
        if not self.tt:
            self.tt=0

        self.TTBE.append(self.tt*self.BULK_TRANS[-1])
        self.TT = (self.TT * self.tt*self.BULK_TRANS[-1])

        return None

    def RestoreData(self):
        """RestoreData.
        """
        for ii in range(0, self.n):
            self.SDT[ii].RestoreSetup()
        self.SetData()

    def StoreData(self):
        """StoreData.
        """
        for ii in range(0, self.n):
            self.SDT[ii].SaveSetup()
        self.SetData()

    def SetData(self):
        """SetData.
        """
        self.SuTo = SUT(self.SDT)
        self.Object_Num = np.arange(0, self.n, 1)
        self.__SurFuncSuscrip()
        self.Pr3D.Prerequisites3SMath()
        self.INORM = InterNormalCalc(self.SDT, self.TypeTotal, self.Pr3D, self.HS)

    def SetSolid(self):
        """SetSolid.
        """
        self.__SurFuncSuscrip()
        self.Pr3D.Prerequisites3SMath()
        self.Pr3D.Prerequisites3DSolids()
        self.INORM = InterNormalCalc(self.SDT, self.TypeTotal, self.Pr3D, self.HS)
        self.AAA = self.Pr3D.AAA
        self.BBB = self.Pr3D.BBB
        self.DDD = self.Pr3D.DDD
        self.EEE = self.Pr3D.EEE
        self.TRANS_1A = self.Pr3D.TRANS_1A
        self.TRANS_2A = self.Pr3D.TRANS_2A

    def Parax(self, W):
        """Parax.

        Parameters
        ----------
        W :
            W
        """
        N_P = []
        for i in range(0, self.n):
            GlGl = self.GlobGlass[i]
            (NP, AP) = n_wave_dispersion(self.SETUP, GlGl, W)
            N_P.append(NP)
        Prx = ParaxCalc(N_P, self.SDT, self.SuTo, self.n, self.Glass)
        (self.SistemMatrix, self.S_Matrix, self.N_Matrix, self.a, self.b, self.c, self.d, self.EFFL, self.PPA, self.PPP, self.c_p, self.n_p, self.d_p) = Prx
        return Prx

    def TargSurf(self, tgsfP1):
        """TargSurf.

        Parameters
        ----------
        tgsfP1 :
            tgsfP1
        """
        tgsf = (tgsfP1 + 1)
        if (self.n >= tgsf > 0):
            self.Targ_Surf = tgsf
        else:
            self.Targ_Surf = self.n

    def SurfFlat(self, fltsf, Prep=0):
        """SurfFlat.

        Parameters
        ----------
        fltsf :
            fltsf
        Prep :
            Prep
        """
        if (self.n >= fltsf > 0):
            self.SuTo.Surface_Flattener = fltsf
        else:
            self.SuTo.Surface_Flattener = 0
        if (Prep != 0):
            self.Pr3D.Prerequisites3DSolids()

    def IgnoreVignetting(self, Prep=0):
        """IgnoreVignetting.

        Parameters
        ----------
        Prep :
            Prep
        """
        self.INORM.Disable_Inner = 0
        self.Pr3D.Disable_Inner = 0
        self.INORM.ExtraDiameter = 1
        self.Pr3D.ExtraDiameter = 1
        if (Prep != 0):
            self.Pr3D.Prerequisites3DSolids()

    def Vignetting(self, Prep=0):
        """Vignetting.

        Parameters
        ----------
        Prep :
            Prep
        """
        self.INORM.Disable_Inner = 1
        self.Pr3D.Disable_Inner = 1
        self.INORM.ExtraDiameter = 0
        self.Pr3D.ExtraDiameter = 0
        if (Prep != 0):
            self.Pr3D.Prerequisites3DSolids()


    def Trace(self, pS, dC, WaveLength):
        """Trace.

        Parameters
        ----------
        pS :
            pS
        dC :
            dC
        WaveLength :
            WaveLength
        """
        self.__CollectDataInit()
        ResVec = np.asarray(dC)
        RayOrig = np.asarray(pS)
        self.RAY = []
        self.Wave = WaveLength
        self.RAY.append(RayOrig)
        self.__WavePrecalc()
        j = 0
        Glass = self.GlobGlass[j]
        (PrevN, alpha) = (self.N_Prec[j], self.AlphaPrecal[j])
        j = 1
        SIGN = np.ones_like(ResVec)
        while True:
            if (j == self.Targ_Surf):
                break
            j_gg = j
            Glass = self.GlobGlass[j_gg]

            if ((self.Glass[j] != 'NULL') and (self.Glass[j] != 'ABSORB')):
                Proto_pTarget = (np.asarray(RayOrig) + ((np.asarray(ResVec) * 999999999.9) * SIGN))


                # start = timeit.default_timer()

                Output = self.INORM.InterNormal(RayOrig, Proto_pTarget, j, j)

                # stop = timeit.default_timer()
                # execution_time = stop - start
                # self.ExectTime.append(execution_time)








                (SurfHit, SurfNorm, pTarget, GooveVect, HitObjSpace, j) = Output

                if (SurfHit == 0):
                    break

                ImpVec = np.asarray(ResVec)
                (CurrN, alpha) = (self.N_Prec[j_gg], self.AlphaPrecal[j_gg])
                S = ImpVec
                R = SurfNorm
                N = PrevN
                Np = CurrN
                D = GooveVect
                Ord = self.SDT[j].Diff_Ord
                GrSpa = self.SDT[j].Grating_D
                Secuent = 0




                # start = timeit.default_timer()

                (ResVec, CurrN, sign) = self.SDT[j].PHYSICS.calculate(S, R, N, Np, D, Ord, GrSpa, self.Wave, Secuent)

                # stop = timeit.default_timer()
                # execution_time = stop - start
                # self.ExectTime.append(execution_time)




                SIGN = (SIGN * sign)
                Name = self.SDT[j].Name
                RayTraceType = 0
                ValToSav = [Glass, alpha, RayOrig, pTarget, HitObjSpace, SurfNorm, ImpVec, ResVec, PrevN, CurrN, WaveLength, D, Ord, GrSpa, Name, j, RayTraceType]
                self.__CollectData(ValToSav)
                PrevN = CurrN
                RayOrig = pTarget
                self.RAY.append(RayOrig)









            j = (j + 1)

        if (len(self.GLASS) == 0):
            self.__CollectDataInit()
            self.val = 0
            self.__EmptyCollect(RayOrig, ResVec, WaveLength, j)

        self.ray_SurfHits = np.asarray(self.RAY)

        AT = np.transpose(self.ray_SurfHits)
        self.Hit_x = AT[0]
        self.Hit_y = AT[1]
        self.Hit_z = AT[2]



        # execution_time=np.mean(np.asarray(self.ExectTime))
        # print(str(execution_time)," Segmentos")


        self.ExectTime=[]

    def NsTrace(self, pS, dC, WaveLength):
        """NsTrace.

        Parameters
        ----------
        pS :
            pS
        dC :
            dC
        WaveLength :
            WaveLength
        """
        global j_gg
        pS = np.asarray(pS)
        count = 0
        a = 1
        b = 2
        self.__CollectDataInit()
        ResVec = dC
        RayOrig = pS
        self.RAY = []
        self.Wave = WaveLength
        self.RAY.append(RayOrig)
        self.__WavePrecalc()
        j = 0
        Glass = self.GlobGlass[j]
        (PrevN, alpha) = (self.N_Prec[j], self.AlphaPrecal[j])
        j = 0
        SIGN = 1
        while True:
            if (j == self.Targ_Surf):
                break
            (a, b, c, PreSurfHit) = self.__NonSequentialChooser(SIGN, RayOrig, ResVec, j)
            if (PreSurfHit == 0):
                break
            if (a < b):
                j_gg = b
            if (a > b):
                j_gg = (b - 1)
                if (self.Glass[(b - 1)] == 'MIRROR'):
                    j_gg = b
                if (self.Glass[b] == 'MIRROR'):
                    j_gg = b
            if (a == b):
                j_gg = (b - 1)
            j = b
            jj = b
            j = self.GlassOnSide[j]
            j_gg = self.GlassOnSide[j_gg]
            Glass = self.GlobGlass[j_gg]
            if ((j == 0) or (count > 20) or (a == self.n)):
                break
            if ((self.Glass[j] != 'NULL') and (self.Glass[j] != 'ABSORB')):
                Proto_pTarget = (np.asarray(RayOrig) + ((np.asarray(ResVec) * 999999999.9) * SIGN))
                Output = self.INORM.InterNormal(RayOrig, Proto_pTarget, j, jj)
                (SurfHit, SurfNorm, pTarget, GooveVect, HitObjSpace, j) = Output
                if (SurfHit == 0):
                    break
                ImpVec = np.asarray(ResVec)
                (CurrN, alpha) = (self.N_Prec[j_gg], self.AlphaPrecal[j_gg])
                S = np.asarray(ImpVec)
                R = np.asarray(SurfNorm)
                N = PrevN
                Np = CurrN
                if ((self.SDT[j].Solid_3d_stl == 'None') and (self.TypeTotal[jj] == 1)):
                    if (N == 1):
                        Np = CurrN
                    else:
                        N = CurrN
                        Np = 1
                D = GooveVect
                Ord = self.SDT[j].Diff_Ord
                GrSpa = self.SDT[j].Grating_D
                Secuent = 0
                (ResVec, CurrN, sign) = self.SDT[j].PHYSICS.calculate(ResVec, R, N, Np, D, Ord, GrSpa, self.Wave, Secuent)
                (Rp0, Rs0, Tp0, Ts0) = FresnelEnergy(self.Glass[j], N, Np, ResVec, R, ResVec, self.SETUP, self.Wave)
                self.tt = 1.0
                if (self.Glass[j] != 'MIRROR'):
                    self.tt = ((Tp0 + Ts0) / 2.0)
                PROB = prob(self.tt)[0]
                if (PROB > 0):
                    Secuent = 1
                    (ResVec, CurrN, sign) = self.SDT[j].PHYSICS.calculate(ResVec, R, N, Np, D, Ord, GrSpa, self.Wave, Secuent)
                SIGN = (SIGN * sign)
                Name = self.SDT[j].Name
                RayTraceType = 1
                ValToSav = [Glass, alpha, RayOrig, pTarget, HitObjSpace, SurfNorm, ImpVec, ResVec, PrevN, CurrN, WaveLength, D, Ord, GrSpa, Name, j, RayTraceType]
                self.__CollectData(ValToSav)
                if (a == b):
                    PrevN = PrevN
                else:
                    PrevN = CurrN
                RayOrig = pTarget
                self.RAY.append(RayOrig)
            count = (count + 1)
        if (len(self.GLASS) == 0):
            self.__CollectDataInit()
            self.val = 0
            self.__EmptyCollect(pS, dC, WaveLength, j)
        self.ray_SurfHits = np.asarray(self.RAY)
        AT = np.transpose(self.ray_SurfHits)
        self.Hit_x = AT[0]
        self.Hit_y = AT[1]
        self.Hit_z = AT[2]

    def FastTrace(self, pS, dC, WaveLength):
        """Trace.

        Parameters
        ----------
        pS :
            pS
        dC :
            dC
        WaveLength :
            WaveLength
        """
        self.__CollectDataInit()
        ResVec = np.asarray(dC)
        RayOrig = np.asarray(pS)
        self.RAY = []
        self.Wave = WaveLength
        self.RAY.append(RayOrig)
        self.__WavePrecalc()
        j = 0
        Glass = self.GlobGlass[j]
        PrevN = self.N_Prec[j]
        j = 1
        SIGN = np.ones_like(ResVec)

        while True:
            if (j == self.Targ_Surf):
                break

            j_gg = j
            Glass = self.GlobGlass[j_gg]

            if ((self.Glass[j] != 'NULL') and (self.Glass[j] != 'ABSORB')):
                Proto_pTarget = (np.asarray(RayOrig) + ((ResVec * 999999999.9) * SIGN))

                Output = self.INORM.InterNormalFast(RayOrig, Proto_pTarget, j)

                (SurfHit, SurfNorm, pTarget, GooveVect, HitObjSpace, j) = Output

                if (SurfHit != 0):

                    ImpVec = ResVec
                    CurrN = self.N_Prec[j_gg]
                    S = ImpVec
                    R = SurfNorm
                    N = PrevN
                    Np = CurrN
                    D = GooveVect
                    Ord = self.SDT[j].Diff_Ord
                    GrSpa = self.SDT[j].Grating_D

                    (ResVec, CurrN, sign) = self.SDT[j].PHYSICS.calculate(S, R, N, Np, D, Ord, GrSpa, self.Wave, 0)

                    SIGN = (SIGN * sign)

                    PrevN = CurrN
                    RayOrig = pTarget
                    self.RAY.append(RayOrig)

                else:
                    break

            j = (j + 1)


        self.ray_SurfHits = np.asarray(self.RAY)


        self.ExectTime=[]
