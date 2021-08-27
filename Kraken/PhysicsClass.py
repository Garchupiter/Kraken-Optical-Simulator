# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 17:13:41 2021

@author: JOELHERRERAVAZQUEZ

"""
import numpy as np

from .Physics import *


##############################################################################

# if self.val == 1:
#             Rp, Rs, Tp, Ts = FresnelEnergy(Glass, PrevN, CurrN, ImpVec, SurfNorm, ResVec, self.SETUP, self.Wave)
#         else:
#             Rp, Rs, Tp, Ts = 0, 0, 0, 0

#         self.RP.append(Rp)
#         self.RS.append(Rs)
#         self.TP.append(Tp)
#         self.TS.append(Ts)

#         if Glass == "MIRROR":
#             tt = 1.0 * (Rp + Rs) / 2.0
#         if Glass != "MIRROR":
#             tt = IT * (Tp + Ts) / 2.0


# if vidrio != "MIRROR":
#     Rp, Rs, Tp, Ts = fresnell_dielectric(NP, NC, ImpVec, SurfNorm, ResVec)

# if vidrio == "MIRROR":
#     n_metal = np.interp(Wave, SETUP.W_alum, SETUP.N_alum)
#     k_complex = np.interp(Wave, SETUP.W_alum, SETUP.K_alum)

#     Rp, Rs, Tp, Ts = fresnell_metal(NP, n_metal, k_complex, ImpVec, SurfNorm)

# return Rp, Rs, Tp, Ts


##############################################################################
class snell_refraction_vector_physics:
    def __init__(self):
        """snell_refraction_vector_physics."""

        """__init__."""
        pass

    def calculate(self, s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, Secuen):
        """calculate.

        :param s1_norm:
        :param Nsurf_norm:
        :param n1:
        :param n2:
        :param empty1:
        :param empty2:
        :param empty3:
        :param empty4:
        :param Secuen:
        """

        s1_norm = np.asarray(s1_norm)
        Nsurf_norm = np.asarray(Nsurf_norm)

        cos = np.dot(Nsurf_norm, s1_norm) / (np.linalg.norm(Nsurf_norm) * np.linalg.norm(s1_norm))
        # print(cos, " ............")
        # print(np.arccos(cos))

        if cos < -1.:
            ang = 180
        else:
            ang = np.rad2deg(np.arccos(cos))
        # print(cos, ang)
        if ang <= 90.0:
            Nsurf_norm = -Nsurf_norm

        Nsurf_Cros_s1 = np.cross(Nsurf_norm, s1_norm)
        SIGN = 1

        if n2 == -1:
            n2 = -n1
            SIGN = -1

        NN = n1 / n2
        R = (NN * NN) * np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1)

        if Secuen == 1:
            R = 2.0

        if R > 1:
            # if Secuen == 0:
            #     print(" - Total internal reflection -")


            n2 = -n1
            NN = n1 / n2
            R = (NN * NN) * np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1)
            SIGN = -1

        s2 = NN * (np.cross(Nsurf_norm, np.cross(-Nsurf_norm, s1_norm))) - Nsurf_norm * np.sqrt(
            1.0 - ((NN * NN) * np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1)))

        orizontal = [0, 0, 1]
        cos = np.dot(orizontal, s2) / (np.linalg.norm(orizontal) * np.linalg.norm(s2))
        ang = np.rad2deg(np.arccos(cos))

        # https://www.embibe.com/study/examples-on-vector-form-of-snell-s-law-concept
        # http://www.starkeffects.com/snells-law-vector.shtml#:~:text=Snell's%20Law%20in%20Vector%20Form&text=Since%20the%20incident%20ray%2C%20the,yourself%20to%20fit%20the%20equation.
        # https://arachnoid.com/OpticalRayTracer_technical/resources/opticalraytracer_technical.pdf

        # verificacion de que los tres vectores son coplanares
        # https://testbook.com/blog/coplanar-vectors/
        # print(np.abs((np.dot(s2,np.cross(Nsurf_norm,s1_norm)))))

        return s2, np.abs(n2), SIGN


##############################################################################
class paraxial_exact_physics:
    def __init__(self):
        """paraxial_exact_physics."""

        """__init__."""
        pass

    def calculate(self, s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, empty5):
        """calculate.

        :param s1_norm:
        :param Nsurf_norm:
        :param n1:
        :param n2:
        :param empty1:
        :param empty2:
        :param empty3:
        :param empty4:
        :param empty5:
        """
        return Nsurf_norm, 1, 1


##############################################################################

class diffraction_grating_physics:
    def __init__(self):
        """diffraction_grating_physics."""

        """__init__."""
        pass

    def calculate(self, S, R, N, Np, D, Ord, d, W, Secuent):
        """calculate.

        :param S:
        :param R:
        :param N:
        :param Np:
        :param D:
        :param Ord:
        :param d:
        :param W:
        :param Secuent:
        """
        # lamb=self.Wave
        lamb = W
        RefRef = -1 * np.sign(Np)
        SIGN = 1
        print(RefRef, Np)

        if Np == -1:
            Np = np.abs(N)

        # RefRef=SIGN

        # U. W. Ludwig, "Generalized grating ray-tracing equations," J. Opt. Soc. Am. 63, 1105-1107 (1973)
        # S= Incident light vector
        # R= Normal to surface
        # D= normal to ruling
        # RefRef =1 for reflection =-1 for refraction
        # d= lines separation

        mu = N / Np
        T = Ord * lamb / (Np * d)

        V = mu * np.dot(S, R)  # (i_r*S[0]+j_r*S[1]+k_r*S[2])
        W = (mu * mu) - 1 + (T * T) - (2.0 * mu * T) * np.dot(S, D)  # (i_d*S[0]+j_d*S[1]+k_d*S[2])

        # Q=np.sqrt(V*V-W)-V  #Reflection
        # Q=-np.sqrt(V*V-W)-V   #refraction
        Q = RefRef * np.sqrt(V * V - W) - V

        # E=Q*Q+2*V*Q+W
        S = np.asarray(S)
        Sp = (mu * S) - (T * D) + (Q * R)
        Sp = Sp / np.linalg.norm(Sp)
        return Sp, np.abs(Np), SIGN

##############################################################
