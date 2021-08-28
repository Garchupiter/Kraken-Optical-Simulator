
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
        """
        s1_norm = np.asarray(s1_norm)
        Nsurf_norm = np.asarray(Nsurf_norm)
        cos = (np.dot(Nsurf_norm, s1_norm) / (np.linalg.norm(Nsurf_norm) * np.linalg.norm(s1_norm)))
        if (cos < (- 1.0)):
            ang = 180
        else:
            ang = np.rad2deg(np.arccos(cos))
        if (ang <= 90.0):
            Nsurf_norm = (- Nsurf_norm)
        Nsurf_Cros_s1 = np.cross(Nsurf_norm, s1_norm)
        SIGN = 1
        if (n2 == (- 1)):
            n2 = (- n1)
            SIGN = (- 1)
        NN = (n1 / n2)
        R = ((NN * NN) * np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1))
        if (Secuen == 1):
            R = 2.0
        if (R > 1):
            n2 = (- n1)
            NN = (n1 / n2)
            R = ((NN * NN) * np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1))
            SIGN = (- 1)
        s2 = ((NN * np.cross(Nsurf_norm, np.cross((- Nsurf_norm), s1_norm))) - (Nsurf_norm * np.sqrt((1.0 - ((NN * NN) * np.dot(Nsurf_Cros_s1, Nsurf_Cros_s1))))))
        orizontal = [0, 0, 1]
        cos = (np.dot(orizontal, s2) / (np.linalg.norm(orizontal) * np.linalg.norm(s2)))
        ang = np.rad2deg(np.arccos(cos))
        return (s2, np.abs(n2), SIGN)

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
        RefRef = ((- 1) * np.sign(Np))
        SIGN = 1
        print(RefRef, Np)
        if (Np == (- 1)):
            Np = np.abs(N)
        mu = (N / Np)
        T = ((Ord * lamb) / (Np * d))
        V = (mu * np.dot(S, R))
        W = ((((mu * mu) - 1) + (T * T)) - (((2.0 * mu) * T) * np.dot(S, D)))
        Q = ((RefRef * np.sqrt(((V * V) - W))) - V)
        S = np.asarray(S)
        Sp = (((mu * S) - (T * D)) + (Q * R))
        Sp = (Sp / np.linalg.norm(Sp))
        return (Sp, np.abs(Np), SIGN)

