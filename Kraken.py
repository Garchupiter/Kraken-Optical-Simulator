"""
Notas realizadas:

"""

import numpy as np

import os
import sys
import random
#import library
currentDirectory = os.getcwd()
sys.path.insert(1, currentDirectory + '/library')

from physics import FresnelEnergy, fresnell_dielectric, fresnell_metal, n_wave_dispersion, ParaxCalc
from display import display3d, display2d
from physics_class import snell_refraction_vector_physics, paraxial_exact_physics, diffraction_grating_physics
from Surf_class import surf
from pupil_tool import pupilcalc, SolveVectCross
from raykeeper import raykeeper
from System_tools import load_alluminum_complex, load_Catalog
from Kraken_setup_class import *
from Surf_tools import surface_tools as SUT
from prerequisites3D import Prerequisites
from HitOnSurf import Hit_Solver
from InterNormalCalc import InterNormalCalc
from WavefrontFit import Zernike_Fitting
from SeidelTool import Seidel

rute = currentDirectory


##############################################################################


def prov(pro):
    a_list = [0, 1]
    prob = pro
    distribution = [prob, 1.0 - prob]
    random_number = random.choices(a_list, distribution)
    return random_number


###########################################################################

class system:

    def __init__(self, SurfData, KN_Setup):

        self.SDT = SurfData

        self.update = False
        self.S_Matrix = []
        self.N_Matrix = []
        self.SistemMatrix = np.matrix([[1.0, 0.0], [0.0, 1.0]])
        self.a, self.b, self.c, self.d, self.EFFL, self.PPA, self.PPP = 0., 0., 0., 0., 0., 0., 0.

        self.SETUP = KN_Setup

        self.n = len(self.SDT)
        self.SuTo = SUT(self.SDT)

        self.Object_Num = np.arange(0, self.n, 1)

        # self.System_stat = self.__PrereqStatus()

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

        self.PreWave = -0.000000001  # Wave lenght for memory

        self.AAA = self.Pr3D.AAA
        self.BBB = self.Pr3D.BBB
        self.DDD = self.Pr3D.DDD
        self.EEE = self.Pr3D.EEE
        self.GlassOnSide = self.Pr3D.GlassOnSide
        self.side_number = self.Pr3D.side_number
        self.TypeTotal = self.Pr3D.TypeTotal
        self.TRANS_1A = self.Pr3D.TRANS_1A

        # print(self.TRANS_1A[2])
        self.TRANS_2A = self.Pr3D.TRANS_2A

        self.HS = Hit_Solver(self.SDT)
        self.INORM = InterNormalCalc(self.SDT, self.TypeTotal, self.Pr3D, self.HS)
        self.INORM.Disable_Inner = 1
        self.Pr3D.Disable_Inner = 1

        self.c_p, self.n_p, self.d_p = 0, 0, 0

    def __SurFuncSuscrip(self):
        for i in range(0, self.n):
            self.SDT[i].build_surface_function()
            # print(i, " in __SurFuncSuscrip")

    def __PrerequisitesGlass(self):
        self.Glass = []
        self.GlobGlass = []

        for i in range(0, self.n):
            if type(self.SDT[i].Glass) == type(1.0):
                self.SDT[i].Glass = "AIR_" + str(self.SDT[i].Glass)
                print(self.SDT[i].Glass)
            self.Glass.append(self.SDT[i].Glass.replace(" ", ""))
            self.GlobGlass.append(self.SDT[i].Glass.replace(" ", ""))

            if self.GlobGlass[i] == "NULL":
                self.GlobGlass[i] = self.GlobGlass[i - 1]

    ###########################################################################

    def __NonSequentialChooserToot(self, A_RayOrig, A_Proto_pTarget, k):
        ng = self.GlassOnSide[k]
        A_Glass = self.SDT[ng].Glass
        if A_Glass == "NULL":
            distance = 99999999999999.9
        else:

            A_pTarget, A_SurfHit = self.EEE[k].ray_trace(A_RayOrig, A_Proto_pTarget)

            if len(A_SurfHit) == 0:
                distance = 99999999999999.9
            else:
                s = 0
                h = []
                for f in A_SurfHit:
                    PD = np.asarray(A_pTarget[s]) - np.asarray(A_RayOrig)
                    distance = np.linalg.norm(PD)
                    if np.abs(distance) < 0.05:
                        distance = 99999999999999.9
                    h.append(distance)
                    s = s + 1
                distance = np.min(np.asarray(h))
        return distance

    ###########################################################################
    def __NonSequentialChooser(self, SIGN, A_RayOrig, ResVec, j):

        chooser = []
        [SLL, SMM, SNN] = ResVec

        distance = 99999999999999.9
        # ZZ=distance*SIGN
        [Px, Py, Pz] = A_RayOrig
        [LL, MM, NN] = ResVec

        A_Proto_pTarget = np.asarray(A_RayOrig) + (np.asarray(ResVec) * 999999999.9 * SIGN)

        for k in range(1, len(self.EEE)):
            distance = self.__NonSequentialChooserToot(A_RayOrig, A_Proto_pTarget, k)

            chooser.append(distance)

        chooser = np.asarray(chooser)
        # print("asdasdasdadsads",chooser)
        jj = np.argmin(chooser) + 1
        chooser[jj - 1] = 99999999999999.9
        kk = np.argmin(chooser) + 1

        A_pTarget, A_SurfHit = self.EEE[int(jj)].ray_trace(A_RayOrig, A_Proto_pTarget)

        PRR = np.shape(A_SurfHit)[0]

        return int(j), int(jj), int(kk), PRR

    ########################################################################

    def __CollectDataInit(self):
        self.val = 1
        self.SURFACE = []
        self.NAME = []
        self.GLASS = []
        self.S_XYZ = []

        self.T_XYZ = []
        self.XYZ = []
        self.XYZ.append([0, 0, 0])

        self.OST_XYZ = []

        self.S_LMN = []  # LMN of surface normal

        self.LMN = []  # LMN of incident ray

        self.R_LMN = []  # LMN of resultant ray

        self.N0 = []
        self.N1 = []
        self.WAV = 1.0

        self.G_LMN = []  # Grating lines direction

        self.ORDER = []
        self.GRATING = []
        self.DISTANCE = []
        self.OP = []
        self.TOP_S = []  # Total optical path step
        self.TOP = 0
        self.ALPHA = []
        self.BULK_TRANS = []
        self.RP = []
        self.RS = []
        self.TP = []
        self.TS = []
        self.TTBE = []  # Total transmision by element
        self.TT = 1.0  # total transmission acumulative

        return None

    #############################################################################

    def __EmptyCollect(self, pS, dC, WaveLength, j):
        Empty = np.asarray([])
        ValToSav = [Empty,
                    Empty,
                    pS,
                    pS,
                    Empty,
                    Empty,
                    dC,
                    Empty,
                    Empty,
                    Empty,
                    WaveLength,
                    Empty,
                    Empty,
                    Empty,
                    Empty,
                    j]
        self.__CollectData(ValToSav)

    def __WavePrecalc(self):
        if self.Wave != self.PreWave:  # Si la longitud de onda es diferente a la previa calcula indices
            self.N_Prec = []  # Indice de refracción precalculado
            self.AlphaPrecal = []  # Alpha o absorción precalculada
            self.PreWave = self.Wave

            for i in range(0, self.n):
                NP, AP = n_wave_dispersion(self.SETUP, self.GlobGlass[i], self.Wave)
                self.N_Prec.append(NP)
                self.AlphaPrecal.append(AP)

            # self.Parax(self.Wave)
            # print(self.Wave)
            # Prx=ParaxCalc(self.N_Prec, self.SDT, self.SuTo, self.n, self.Glass)
            # self.SistemMatrix, self.S_Matrix, self.N_Matrix, self.a, self.b, self.c, self.d, self.EFFL, self.PPA, self.PPP=Prx

    def __CollectData(self, ValToSav):
        [Glass,
         alpha,
         RayOrig,
         pTarget,
         HitObjSpace,
         SurfNorm,
         ImpVec,
         ResVec,
         PrevN,
         CurrN,
         WaveLength,
         D,
         Ord,
         GrSpa,
         Name,
         j] = ValToSav

        global tt
        self.SURFACE.append(j)
        self.NAME.append(Name)
        self.GLASS.append(Glass)

        self.S_XYZ.append(RayOrig)
        self.T_XYZ.append(pTarget)
        self.XYZ[0] = self.S_XYZ[0]
        self.XYZ.append(pTarget)

        self.OST_XYZ.append(HitObjSpace)

        p = np.asarray(RayOrig) - np.asarray(pTarget)
        dist = np.linalg.norm(p)
        self.DISTANCE.append(dist)
        self.OP.append(dist * PrevN)
        self.TOP = self.TOP + (dist * PrevN)
        self.TOP_S.append(self.TOP)
        self.ALPHA.append(alpha)
        IT = np.exp(-alpha * dist)
        self.BULK_TRANS.append(IT)

        self.S_LMN.append(SurfNorm)
        self.LMN.append(ImpVec)
        self.R_LMN.append(ResVec)

        self.N0.append(PrevN)
        self.N1.append(CurrN)
        self.WAV = WaveLength

        self.G_LMN.append(D)

        self.ORDER.append(Ord)
        self.GRATING.append(GrSpa)
        if self.val == 1:
            Rp, Rs, Tp, Ts = FresnelEnergy(Glass, PrevN, CurrN, ImpVec, SurfNorm, ResVec, self.SETUP, self.Wave)
        else:
            Rp, Rs, Tp, Ts = 0, 0, 0, 0

        self.RP.append(Rp)
        self.RS.append(Rs)
        self.TP.append(Tp)
        self.TS.append(Ts)

        if Glass == "MIRROR":
            tt = 1.0 * (Rp + Rs) / 2.0
        if Glass != "MIRROR":
            tt = IT * (Tp + Ts) / 2.0

        self.TTBE.append(tt)
        self.TT = self.TT * tt
        return None

    def ResetData(self):
        self.SuTo = SUT(self.SDT)
        self.Object_Num = np.arange(0, self.n, 1)
        self.__SurFuncSuscrip()
        # self.Pr3D = Prerequisites(self.SDT, self.SuTo)
        # self.__PrerequisitesGlass()
        self.Pr3D.Prerequisites3SMath()

        self.INORM = InterNormalCalc(self.SDT, self.TypeTotal, self.Pr3D, self.HS)

    def ResetSolid(self):
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

        ##  "sest" mean system estatus
        # sest, sest_validation = self.ModificationDetector()
        # if sest_validation==True:
        #     self.System_stat = self.System_stat_at_this_time
        #     self.__SecuentialStatusChange()
        #     self.__SurFuncSuscrip()
        # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--

        N_P = []  # Indice de refracción precalculado

        for i in range(0, self.n):
            GlGl = self.GlobGlass[i]
            NP, AP = n_wave_dispersion(self.SETUP, GlGl, W)
            N_P.append(NP)

        Prx = ParaxCalc(N_P, self.SDT, self.SuTo, self.n, self.Glass)
        # SistemMatrix, S_Matrix, N_Matrix, a, b, c, d, EFFL, PPA, PPP=Prx
        self.SistemMatrix, self.S_Matrix, self.N_Matrix, self.a, self.b, self.c, self.d, self.EFFL, self.PPA, self.PPP, self.c_p, self.n_p, self.d_p = Prx
        return Prx

    def TargSurf(self, tgsfP1):

        tgsf = tgsfP1 + 1
        if self.n >= tgsf > 0:
            self.Targ_Surf = tgsf
            # print("End (Target) surface changed")
        else:
            # print("Error: target number is greater than the number of elements, default value is used ")
            self.Targ_Surf = self.n

    def SurfFlat(self, fltsf, Prep=0):

        if self.n >= fltsf > 0:
            self.SuTo.Surface_Flattener = fltsf

            # print("Flat surface defined")
        else:
            # print("Warning: surface number is greater than the number of elements, default value is used ")
            # print("Don't worry if it was set on purpose with -1 to return to default value")
            self.SuTo.Surface_Flattener = 0
        if Prep != 0:
            self.Pr3D.Prerequisites3DSolids()

    def IgnoreVignetting(self, Prep=0):
        self.INORM.Disable_Inner = 0
        self.Pr3D.Disable_Inner = 0

        self.INORM.ExtraDiameter = 1
        self.Pr3D.ExtraDiameter = 1

        if Prep != 0:
            self.Pr3D.Prerequisites3DSolids()

    def Vignetting(self, Prep=0):
        self.INORM.Disable_Inner = 1
        self.Pr3D.Disable_Inner = 1

        self.INORM.ExtraDiameter = 0
        self.Pr3D.ExtraDiameter = 0

        if Prep != 0:
            self.Pr3D.Prerequisites3DSolids()

    def Trace(self, pS, dC, WaveLength):

        self.__CollectDataInit()
        ResVec = np.asarray(dC)
        RayOrig = np.asarray(pS)
        self.RAY = []

        self.Wave = WaveLength
        self.RAY.append(RayOrig)
        self.__WavePrecalc()

        j = 0
        Glass = self.GlobGlass[j]
        PrevN, alpha = self.N_Prec[j], self.AlphaPrecal[j]

        j = 1
        SIGN = np.ones_like(ResVec)

        while True:
            if j == self.Targ_Surf:
                break

            j_gg = j
            Glass = self.GlobGlass[j_gg]

            if self.Glass[j] != "NULL" and self.Glass[j] != "ABSORB":  # Si es diferente de nulo

                Proto_pTarget = np.asarray(RayOrig) + (np.asarray(ResVec) * 999999999.9 * SIGN)

                Output = self.INORM.InterNormal(RayOrig, Proto_pTarget, j, j)
                SurfHit, SurfNorm, pTarget, GooveVect, HitObjSpace, j = Output

                if SurfHit == 0:
                    break

                ImpVec = np.asarray(ResVec)
                CurrN, alpha = self.N_Prec[j_gg], self.AlphaPrecal[j_gg]

                S = ImpVec
                R = SurfNorm
                N = PrevN
                Np = CurrN
                D = GooveVect
                Ord = self.SDT[j].Diff_Ord
                GrSpa = self.SDT[j].Grating_D
                Secuent = 0
                ResVec, CurrN, sign = self.SDT[j].PHYSICS.calculate(S, R, N, Np, D, Ord, GrSpa, self.Wave, Secuent)

                SIGN = SIGN * sign
                Name = self.SDT[j].Name
                ValToSav = [Glass,
                            alpha,
                            RayOrig,
                            pTarget,
                            HitObjSpace,
                            SurfNorm,
                            ImpVec,
                            ResVec,
                            PrevN,
                            CurrN,
                            WaveLength,
                            D,
                            Ord,
                            GrSpa,
                            Name,
                            j]
                self.__CollectData(ValToSav)
                PrevN = CurrN
                RayOrig = pTarget
                self.RAY.append(RayOrig)
            j = j + 1
        ##----------------------------------------------------------------------
        if len(self.GLASS) == 0:
            self.__CollectDataInit()
            self.val = 0

            self.__EmptyCollect(RayOrig, ResVec, WaveLength, j)

        self.ray_SurfHits = np.asarray(self.RAY)
        AT = np.transpose(self.ray_SurfHits)
        self.Hit_x = AT[0]
        self.Hit_y = AT[1]
        self.Hit_z = AT[2]

    ##############################################################################

    def NsTrace(self, pS, dC, WaveLength):
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
        PrevN, alpha = self.N_Prec[j], self.AlphaPrecal[j]

        j = 0
        SIGN = 1

        while True:
            if j == self.Targ_Surf:
                break

            ##--------------------------------------------------

            a, b, c, PreSurfHit = self.__NonSequentialChooser(SIGN, RayOrig, ResVec,
                                                              j)  # Investiga cual de los elementos es tocado por el rayo

            # if b > (self.n-1) or PreSurfHit==0:
            if PreSurfHit == 0:
                break
            if a < b:
                j_gg = b

            if a > b:
                j_gg = b - 1

                if self.Glass[b - 1] == "MIRROR":
                    j_gg = b

                if self.Glass[b] == "MIRROR":
                    j_gg = b

            if a == b:
                j_gg = b - 1

            j = b
            jj = b
            j = self.GlassOnSide[j]

            j_gg = self.GlassOnSide[j_gg]
            Glass = self.GlobGlass[j_gg]

            if j == 0 or count > 20 or a == self.n:
                break

            if self.Glass[j] != "NULL" and self.Glass[j] != "ABSORB":  # Si es diferente de nulo
                Proto_pTarget = np.asarray(RayOrig) + (np.asarray(ResVec) * 999999999.9 * SIGN)
                Output = self.INORM.InterNormal(RayOrig, Proto_pTarget, j, jj)
                SurfHit, SurfNorm, pTarget, GooveVect, HitObjSpace, j = Output

                if SurfHit == 0:
                    break

                ImpVec = np.asarray(ResVec)
                CurrN, alpha = self.N_Prec[j_gg], self.AlphaPrecal[j_gg]
                S = np.asarray(ImpVec)
                R = np.asarray(SurfNorm)
                N = PrevN
                Np = CurrN
                if self.SDT[j].Solid_3d_stl == "None" and self.TypeTotal[jj] == 1:
                    if N == 1:
                        Np = CurrN
                    else:
                        N = CurrN
                        Np = 1

                D = GooveVect
                Ord = self.SDT[j].Diff_Ord
                GrSpa = self.SDT[j].Grating_D
                Secuent = 0
                ResVec, CurrN, sign = self.SDT[j].PHYSICS.calculate(ResVec, R, N, Np, D, Ord, GrSpa, self.Wave, Secuent)
                Rp0, Rs0, Tp0, Ts0 = FresnelEnergy(self.Glass[j], N, Np, ResVec, R, ResVec, self.SETUP, self.Wave)
                tt = 1.0
                if self.Glass[j] != "MIRROR":
                    tt = (Tp0 + Ts0) / 2.0

                PROB = prov(tt)[0]

                # print(self.Glass[j], tt, PROB)

                if PROB > 0:
                    Secuent = 1
                    ResVec, CurrN, sign = self.SDT[j].PHYSICS.calculate(ResVec, R, N, Np, D, Ord, GrSpa, self.Wave,
                                                                        Secuent)

                SIGN = SIGN * sign
                Name = self.SDT[j].Name
                ValToSav = [Glass,
                            alpha,
                            RayOrig,
                            pTarget,
                            HitObjSpace,
                            SurfNorm,
                            ImpVec,
                            ResVec,
                            PrevN,
                            CurrN,
                            WaveLength,
                            D,
                            Ord,
                            GrSpa,
                            Name,
                            j]
                self.__CollectData(ValToSav)

                if a == b:
                    PrevN = PrevN
                else:
                    PrevN = CurrN
                RayOrig = pTarget
                self.RAY.append(RayOrig)

            count = count + 1

        ##----------------------------------------------------------------------

        if len(self.GLASS) == 0:
            self.__CollectDataInit()
            self.val = 0

            self.__EmptyCollect(pS, dC, WaveLength, j)

        self.ray_SurfHits = np.asarray(self.RAY)
        AT = np.transpose(self.ray_SurfHits)
        self.Hit_x = AT[0]
        self.Hit_y = AT[1]
        self.Hit_z = AT[2]
