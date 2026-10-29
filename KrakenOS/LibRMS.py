
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
    # The RMS radius has no zero for a real bundle: asking fsolve for one returned wherever it
    # stopped. Moved by delta_Z a ray
    # lands at x + delta_Z * L / N, so the squared RMS radius is a parabola in delta_Z and its
    # minimum is exact. Parallel rays have no best focus: 0.
    X, Y, L, M, N = (np.asarray(value, dtype=float).reshape(-1) for value in (X, Y, L, M, N))
    u, v = L / N, M / N
    x1, y1 = X - np.mean(X), Y - np.mean(Y)
    u1, v1 = u - np.mean(u), v - np.mean(v)
    spread = float(np.mean(u1 * u1 + v1 * v1))
    if not np.isfinite(spread) or spread <= 0.0:
        return 0.0
    return float(-np.mean(x1 * u1 + y1 * v1) / spread)

