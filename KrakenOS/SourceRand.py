
from scipy.stats import uniform
import random
import numpy as np

class SourceRnd():
    """SourceRnd.

        SourceRnd.fun = Python function of angle to cover
        SourceRnd.field = Angle to cover from pole (Deg)
        SourceRnd.type = 0 for circle 1 to square
        SourceRnd.dim = Lateral size
        SourceRnd.num = Number of rays

    """


    def __init__(self):
        """__init__.
        """
        self.fun = 0
        self.field = 89.9
        self.type = 0
        self.dim = 10
        self.num = 100

    def rays(self):
        """rays.
        """

        def ff(x):
            """ff.

            Parameters
            ----------
            x :
                x
            """
            return (x / x)
        f = self.fun
        if (f == 0):
            f = ff
        r_source = self.dim
        lim = np.deg2rad(self.field)
        A = np.linspace((lim / self.num), (lim - (lim / self.num)), self.num)
        prob = (f(A) * ((2 * np.pi) * np.sin(A)))
        prob = (prob - np.min(prob))
        prob = (prob / np.max(prob))
        Rand_num = random.choices(A, prob, k=self.num)
        Delta = ((lim / self.num) * (uniform.rvs(size=self.num, loc=0, scale=2) - 1))
        Theta = (Rand_num + Delta)
        Phi = uniform.rvs(size=self.num, loc=0, scale=(2.0 * np.pi))
        L = (np.sin(Theta) * np.cos(Phi))
        M = (np.sin(Theta) * np.sin(Phi))
        N = np.cos(Theta)
        if (self.type == 0):
            AR = np.linspace((r_source / self.num), (r_source - (r_source / self.num)), self.num)
            prob_R = ((2 * np.pi) * AR)
            lt0 = np.argwhere((prob_R < 0))
            prob_R = (prob_R - np.min(prob_R))
            prob_R = (prob_R / np.max(prob_R))
            Rand_num_R = random.choices(AR, prob_R, k=self.num)
            Delta_R = ((lim / self.num) * (uniform.rvs(size=self.num, loc=0, scale=2) - 1))
            R = (Rand_num_R + Delta_R)
            Phi_R = uniform.rvs(size=self.num, loc=0, scale=(2.0 * np.pi))
            X = (R * np.cos(Phi_R))
            Y = (R * np.sin(Phi_R))
        if (self.type == 1):
            X = (r_source * (uniform.rvs(size=self.num, loc=0, scale=2) - 1))
            Y = (r_source * (uniform.rvs(size=self.num, loc=0, scale=2) - 1))
        Z = np.zeros_like(X)
        return (L, M, N, X, Y, Z)

