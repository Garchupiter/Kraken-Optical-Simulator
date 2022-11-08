
import numpy as np
from .Physics import *

class snell_refraction_vector_physics():
    """snell_refraction_vector_physics.
    """


    def __init__(self):
        """__init__.
        """
        pass

    def calculate(self, s1_norm, Nsurf_norm, n1, n2, empty1, empty2, empty3, empty4, Secuen):
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
        else:
            ang = 180 - ang
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
            # Replexion total interna
            n2 = - n1
            NN = n1 / n2
            SIGN = -1.0

        c1 = np.dot(Nv,Iv)
        if c1 < 0.0 :
            c1 = np.dot(-Nv,Iv)
        IP = ((NN**2)*(1-(c1**2.0)))

        c2 = np.sqrt(1.0-IP)
        T = NN * Iv + ((( NN * c1 ) - c2)) * Nv

        return T, np.abs(n2), SIGN, ang

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
        return (Nsurf_norm, 1, 1, 0)

class diffraction_grating_physics():
    """diffraction_grating_physics.
    """

    def __init__(self):
        """__init__.
        """
        pass

    def calculate(self, S, R, N, Np, P, Ord, d, W, Secuent):

        cos=np.dot(R,S)
        Ang=np.rad2deg(np.arccos( cos))

        ang = Ang

        if (cos < (- 1.0)):
            ang = 180.0
        else:
            ang = np.rad2deg(np.arccos(cos))

        if (ang >= 90.0):
            ang = 180 - ang


        if np.abs(Ang) < 90:
            R = -R


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

        D = np.cross(R, P)

        lamb = W
        RefRef = ((- 1.) * np.sign(Np))
        SIGN = 1.

        if (Np == (- 1.)):
            Np = np.abs(N)

        mu = (N / Np)
        T = ((Ord * lamb) / (Np * d))

        V = (mu * np.dot(R,S))
        W = ((((mu ** 2.0) - 1.) + (T ** 2.0)) - (((2.0 * mu) * T) * np.dot(D,S)))

        Q1 = ( np.sqrt((V**2 - W))) - V
        Q2 = (-np.sqrt((V**2 - W))) - V

        if RefRef == 1:
            if Q1 > Q2:
                Q = Q1
            else:
                Q = Q2

        if RefRef == -1:
            if Q1 < Q2:
                Q = Q1
            else:
                Q = Q2

        S = np.asarray(S)
        Sp = (((mu * S) - (T * D)) + (Q * R))
        SIGN = SIGN*-1*RefRef
        return (Sp, np.abs(Np), SIGN, ang)

