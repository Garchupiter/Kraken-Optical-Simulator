
import numpy as np
import scipy

def R_RMS_delta(delta_Z, L, M, N, X, Y):
    """R_RMS_delta.

    Parameters
    ----------
    delta_Z :
        delta_Z
    L :
        L
    M :
        M
    N :
        N
    X :
        X
    Y :
        Y
    """
    X = (((L / N) * delta_Z) + X)
    Y = (((M / N) * delta_Z) + Y)
    cenX = np.mean(X)
    cenY = np.mean(Y)
    x1 = (X - cenX)
    y1 = (Y - cenY)
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    return R_RMS

def RMS(X, Y, Z, L, M, N):
    """RMS.

    Parameters
    ----------
    X :
        X
    Y :
        Y
    Z :
        Z
    L :
        L
    M :
        M
    N :
        N
    """
    cenX = np.mean(X)
    cenY = np.mean(Y)
    x1 = (X - cenX)
    y1 = (Y - cenY)
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    return (R_RMS, cenX, cenY)

def BestFocus(X, Y, Z, L, M, N):
    """BestFocus.

    Parameters
    ----------
    X :
        X
    Y :
        Y
    Z :
        Z
    L :
        L
    M :
        M
    N :
        N
    """
    delta_Z = 0
    ZZ = (L, M, N, X, Y)
    v = scipy.optimize.fsolve(R_RMS_delta, delta_Z, args=ZZ)
    return v[0]

