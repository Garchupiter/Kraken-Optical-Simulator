
import numpy as np
from .Physics import *
from numba import jit

# from numba import double, autojit

class snell_refraction_vector_physics():
    """snell_refraction_vector_physics.
    """


    def __init__(self):
        """__init__.
        """
        pass

    def calculate(self, s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, Secuen):
    #     """calculate.

    #     Parameters
    #     ----------
    #     s1_norm :
    #         s1_norm
    #     Nsurf_norm :
    #         Nsurf_norm
    #     n1 :
    #         n1
    #     n2 :
    #         n2
    #     empty1 :
    #         empty1
    #     empty2 :
    #         empty2
    #     empty3 :
    #         empty3
    #     empty4 :
    #         empty4
    #     Secuen :
    #         Secuen
    #     """

    #     cos = (np.dot(Nsurf_norm, s1_norm))

    #     if (cos < (- 1.0)):
    #         ang = 180
    #     else:
    #         ang = np.rad2deg(np.arccos(cos))

    #     if (ang <= 90.0):
    #         Nsurf_norm = (- Nsurf_norm)

    #     Nsurf_Cros_s1 = np.cross(Nsurf_norm, s1_norm)
    #     SIGN = 1

    #     if (n2 == (- 1)):
    #         n2 = (- n1)
    #         SIGN = (- 1)

    #     NN = (n1 / n2)
    #     d22=np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1)
    #     R = ((NN * NN) * d22)

    #     if (Secuen == 1):
    #         R = 2.0

    #     if (R > 1):
    #         n2 = (- n1)
    #         NN = (n1 / n2)
    #         R = ((NN * NN) * d22)
    #         SIGN = (- 1)

    #     s2 = ((NN * np.cross(Nsurf_norm, np.cross((- Nsurf_norm), s1_norm))) - (Nsurf_norm * np.sqrt((1.0 - ((NN * NN) * d22)))))


    #     cos = np.dot([0, 0, 1], s2)

    #     ang = np.rad2deg(np.arccos(cos))


    #     return (s2, np.abs(n2), SIGN)

        return calculate2(s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, Secuen)


# @jit(nopython=True)
def calculate2(s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, Secuen):
    """calculate.

    Parameters
    ----------
    s1_norm :
        s1_norm
    Nsurf_norm :
        Nsurf_norm
    n1 :
        n1
    n2 :
        n2
    empty1 :
        empty1
    empty2 :
        empty2
    empty3 :
        empty3
    empty4 :
        empty4
    Secuen :
        Secuen
r = rayo
r'=reflejado
S=normal superficie
    """


    Nv = np.asarray(Nsurf_norm)
    Iv = np.asarray(s1_norm)


    """ checking if the normal to the surface is in the
    direction of space where the ray is coming from and
    not where it is going """

    cos=np.dot(Nv, Iv)

    if (cos < (- 1.0)):
        ang = 180.0
    else:
        ang = np.rad2deg(np.arccos(cos))

    if (ang <= 90.0):
        Nv = - Nv

    """----------------------------------------------"""


    Nsurf_Cros_s1 = np.cross(Nv, Iv)
    SIGN = 1.0

    if (n2 == - 1.0):
        n2 = - n1
        SIGN = -1.0

    NN = (n1 / n2)
    d22=np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1)
    R = (NN * NN) * d22

    if (Secuen == 1.0):
        R = 2.0

    if (R > 1.0):
        n2 = - n1
        NN = n1 / n2
        # R = (NN * NN) * d22
        SIGN = - 1.0

    c1 = np.dot(Nv,Iv)
    if c1 < 0.0 :
        c1 = np.dot(-Nv,Iv)
    IP = ((NN**2)*(1-(c1**2.0)))

    c2 = np.sqrt(1.0-IP)
    T = NN * Iv + ((( NN * c1 ) - c2)) * Nv

    return T, np.abs(n2), SIGN


# @jit(nopython=True)
# def calculate2(s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, Secuen):
#     """calculate.

#     Parameters
#     ----------
#     s1_norm :
#         s1_norm
#     Nsurf_norm :
#         Nsurf_norm
#     n1 :
#         n1
#     n2 :
#         n2
#     empty1 :
#         empty1
#     empty2 :
#         empty2
#     empty3 :
#         empty3
#     empty4 :
#         empty4
#     Secuen :
#         Secuen
#     """


#     # Nsurf_norm = 1. * np.asarray(Nsurf_norm1)
#     # s1_norm = 1. * np.asarray(s1_norm1)


#     N1 = [Nsurf_norm[0], Nsurf_norm[1], Nsurf_norm[2]]
#     S1 = [s1_norm[0], s1_norm[1], s1_norm[2]]

#     cos=(N1[0]*S1[0]) + (N1[1]*S1[1]) + (N1[2]*S1[2])

#     # cos = np.dot(N1, S1)

#     if (cos < (- 1.0)):
#         ang = 180.0
#     else:
#         ang = np.rad2deg(np.arccos(cos))

#     if (ang <= 90.0):
#         Nsurf_norm = (- Nsurf_norm)

#     Nsurf_Cros_s1 = np.cross(Nsurf_norm, s1_norm)
#     SIGN = 1.0

#     if (n2 == (- 1.0)):
#         n2 = (- n1)
#         SIGN = (- 1.0)

#     NN = (n1 / n2)
#     d22=np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1)
#     R = ((NN * NN) * d22)

#     if (Secuen == 1.0):
#         R = 2.0

#     if (R > 1.0):
#         n2 = (- n1)
#         NN = (n1 / n2)
#         R = ((NN * NN) * d22)
#         SIGN = (- 1.0)

#     s2 = ((NN * np.cross(Nsurf_norm, np.cross((- Nsurf_norm), s1_norm))) - (Nsurf_norm * np.sqrt((1.0 - ((NN * NN) * d22)))))


#     # cos = np.dot([0.0, 0.0, 1.0], [s2[0], s2[1], s2[2]])
#     # cos = s2[2]



#     # ang = np.rad2deg(np.arccos(cos))


#     return s2, np.abs(n2), SIGN





class paraxial_exact_physics():
    """paraxial_exact_physics.
    """


    def __init__(self):
        """__init__.
        """
        pass

    def calculate(self, s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, empty5):
        """calculate.

        Parameters
        ----------
        s1_norm :
            s1_norm
        Nsurf_norm :
            Nsurf_norm
        n1 :
            n1
        n2 :
            n2
        empty1 :
            empty1
        empty2 :
            empty2
        empty3 :
            empty3
        empty4 :
            empty4
        empty5 :
            empty5
        """
        return (Nsurf_norm, 1, 1)





class diffraction_grating_physics():
    """diffraction_grating_physics.
    """


    def __init__(self):
        """__init__.
        """
        pass

    def calculate(self, S, R, N, Np, D, Ord, d, W, Secuent):
        """calculate.

        Parameters
        ----------
        S :
            S
        R :
            R
        N :
            N
        Np :
            Np
        D :
            D
        Ord :
            Ord
        d :
            d
        W :
            W
        Secuent :
            Secuent
        """
        lamb = W
        RefRef = ((- 1.) * np.sign(Np))
        SIGN = 1.
        print(RefRef, Np)
        if (Np == (- 1.)):
            Np = np.abs(N)
        mu = (N / Np)
        T = ((Ord * lamb) / (Np * d))
        V = (mu * np.dot(S, R))
        W = ((((mu * mu) - 1.) + (T * T)) - (((2.0 * mu) * T) * np.dot(S, D)))
        Q = ((RefRef * np.sqrt(((V * V) - W))) - V)
        S = np.asarray(S)
        Sp = (((mu * S) - (T * D)) + (Q * R))
        Sp = (Sp / np.linalg.norm(Sp))
        return (Sp, np.abs(Np), SIGN)

