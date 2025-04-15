import KrakenOS as Kos
import pkg_resources

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import pickle
import scipy
from scipy.optimize import curve_fit
from pickle import DICT
import copy
###############################################################################

def adjust_tilt(X, Y, Z):
    """
    Fits a plane to the phase map given by arrays X, Y, and Z
    to remove tilt and obtain the adjusted phase map.

    Parameters:
    - X: 2D array with the X coordinates of the phase map.
    - Y: 2D array with the Y coordinates of the phase map.
    - Z: 2D array with phase error values at each point (phase map).

    Returns:
    - Z_adjusted: 2D array with the adjusted phase map, without tilt.
    """

    # Flatten the X, Y, Z arrays to perform the fitting
    x_flat = X.ravel()
    y_flat = Y.ravel()
    z_flat = Z.ravel()

    # Define a plane function (tilt)
    def plane(xy, a, b, c):
        x, y = xy
        return a * x + b * y + c

    # Fit a plane to the Z data as a function of X and Y
    xy = np.vstack((x_flat, y_flat))
    popt, _ = curve_fit(plane, xy, z_flat)

    # Get the fitted plane coefficients
    a, b, c = popt

    # Calculate the fitted plane over the original X, Y coordinates
    plane_adjusted = a * X + b * Y + c

    # Subtract the fitted plane to obtain the phase map without tilt
    Z_adjusted = Z - plane_adjusted

    return Z_adjusted

###############################################################################

def precise_det_exit_pupil(Sys, Pupil, w):
    DeltaField = 0.01
    DFx = [DeltaField, -DeltaField, 0, 0]
    DFy = [0, 0, DeltaField, -DeltaField]

    XA = []
    YA = []
    ZA = []
    LA = []
    MA = []
    NA = []

    for i in range(4):
        Pupil.FieldY = DFx[i]
        Pupil.FieldX = DFy[i]
        Pupil.Pattern()
        XX, YY, ZZ, LL, MM, NN = Pupil.Pattern2FieldPlus()
        Sys.Trace([XX[0], YY[0], ZZ[0]], [LL[0], MM[0], NN[0]], w)
        [Jx, Jy, Jz] = Sys.XYZ[-1]
        [Jl, Jm, Jn] = Sys.LMN[-1]
        XA.append(Jx)
        YA.append(Jy)
        ZA.append(Jz)
        LA.append(Jl)
        MA.append(Jm)
        NA.append(Jn)

    XA = np.asarray(XA)
    YA = np.asarray(YA)
    ZA = np.asarray(ZA)
    LA = np.asarray(LA)
    MA = np.asarray(MA)
    NA = np.asarray(NA)

    POZ = Kos.BestExitPupilPos(XA, YA, ZA, LA, MA, NA)
    return POZ

def R_RMS_delta(Z1, L, M, N, X0, Y0):
    X1 = ((L / N) * Z1) + X0
    Y1 = ((M / N) * Z1) + Y0
    cenX = np.mean(X1)
    cenY = np.mean(Y1)
    x1 = (X1 - cenX)
    y1 = (Y1 - cenY)
    R2 = ((x1 * x1) + (y1 * y1))
    R_RMS = np.sqrt(np.mean(R2))
    return R_RMS

def BestExitPupilPos(X, Y, Z, L, M, N):
    delta_Z = 0
    ZZ = (L, M, N, X, Y)
    vz = scipy.optimize.fsolve(R_RMS_delta, delta_Z, args=ZZ)

    Delta_Z = vz[0]
    X1 = ((L / N) * Delta_Z) + X
    Y1 = ((M / N) * Delta_Z) + Y
    cenX = np.mean(X1)
    cenY = np.mean(Y1)

    return Delta_Z

###############################################################################
def Cat2SurfSimple(cat_dict: DICT):
    '''
    Convert lens catalog dict to surface list

    Parameters
    ----------
    cat_dict : DICT
        Parsed lens catalog dict
    No Desp, Tilts, and no AxisMove; simple surface

    Returns
    -------
    _type_
        surface list
    '''

    surf_name = [surface for surface in cat_dict.keys() if ('SUFR' in surface)]
    surf_list = []

    n = len(surf_name)
    i = 0

    for idx, surface in enumerate(surf_name):
        sf = Kos.surf()

        current_surf = cat_dict[surface]
        sf.Rc = current_surf.get('Rc', 0)
        Thickness = current_surf.get('Thickness', 0)
        if Thickness == 0:
            # if i == 0:
            #     Thickness = 10
            # else:
            Thickness = 0.0001

        sf.Thickness = Thickness
        Diameter = current_surf.get('Diameter', 1)
        if Diameter == 0:
            Diameter = 1
        sf.Diameter = Diameter
        sf.k = current_surf.get('conic', 0)
        sf.AspherData = np.array(current_surf.get('aspherics', [0]*200))
        sf.Glass = current_surf.get('Glass', 'AIR')

        surf_list.append(sf)
        i = i+1

    return surf_list

###############################################################################
def find_intersections(pSource, dCos, x0_v, y0_v, z0_v):
    """
    Finds intersection points between a plane and a series of vectors
    that share the same direction cosines.

    Parameters:
    - pSource: list or array with the origin point of the plane's normal vector [x, y, z].
    - dCos: list or array with the direction cosines of the plane's normal vector [L, M, N].
    - x0_v, y0_v, z0_v: lists or arrays of origin points of other vectors, where each vector
    starts at (x0_v[i], y0_v[i], z0_v[i]).

    Returns:
    - ix, iy, iz: numpy arrays with the x, y, z coordinates of the intersection points.
    """
    # Define the plane with coefficients A, B, C, and D
    x0, y0, z0 = pSource
    L, M, N = dCos
    D = -(L * x0 + M * y0 + N * z0)

    # List to store intersection points
    ix = []
    iy = []
    iz = []

    # Calculate the parameter t for the intersection
    denominator = L * L + M * M + N * N
    for i in range(len(x0_v)):

        t = -(L * x0_v[i] + M * y0_v[i] + N * z0_v[i] + D) / denominator

        # Calculate the intersection point (x, y, z)
        x_inter = x0_v[i] + t * L
        y_inter = y0_v[i] + t * M
        z_inter = z0_v[i] + t * N

        ix.append(x_inter)
        iy.append(y_inter)
        iz.append(z_inter)

    ix = np.asarray(ix)
    iy = np.asarray(iy)
    iz = np.asarray(iz)

    return ix, iy, iz

###############################################################################
def calculate_angles_tx_ty(dCos_new):
    """
    Calculates the Tx and Ty angles that transform the initial vector [0, 0, 1]
    to the vector given by the direction cosines dCos_new = [L, M, N].

    Parameters:
    - dCos_new: list or array with direction cosines [L, M, N].

    Returns:
    - Tx: Rotation angle around the X-axis in radians.
    - Ty: Rotation angle around the Y-axis in radians.
    """
    L, M, N = dCos_new

    # Calculate Ty
    Ty = np.arcsin(L)

    # Calculate Tx using M and N
    Tx = np.arctan2(-M, N)

    return np.rad2deg(Tx), np.rad2deg(Ty)

def RequestRays(Rays, Pupil):
    # Request all ray data on the first surface
    X, Y, Z, L, M, N = Rays.pick(0)
    Lx = (X[3] - X[1]) / 2
    Ly = (Y[4] - Y[2]) / 2
    Cx = X[0]
    Cy = Y[0]

    Pupil.Ptype = "hexapolar"
    Pupil.Pattern()
    Px = (Pupil.Cordx * Lx) + Cx
    Py = (Pupil.Cordy * Ly) + Cy
    Pz = np.zeros_like(Px)

    pSource = [X[0], Y[0], Z[0]]
    dCos = [L[0], M[0], N[0]]

    NPx, NPy, NPz = Kos.find_intersections(pSource, dCos, Px, Py, Pz)
    return NPx, NPy, NPz, L, M, N, Pupil

def CreateRefSphere(Rays, Pupil, SYS, POZ, NPx, NPy, NPz, L, M, N, w, CR_pSource):
    ###############################################################################
    """ For the chief ray in the image plane CRPI """
    RPI_X, RPI_Y, RPI_Z, RPI_L, RPI_M, RPI_N = Rays.pick(-1)

    Xf = RPI_X[0]
    Yf = RPI_Y[0]
    Zf = RPI_Z[0]
    Lf = RPI_L[0]
    Mf = RPI_M[0]
    Nf = RPI_N[0]

    """ Calculate the tilt angles of the reference sphere with the chief ray """
    Tx, Ty = Kos.calculate_angles_tx_ty([Lf, Mf, Nf])

    """ Determine the lateral displacement of the pupil with the intersection
    of the chief ray at POZ """
    DespEsfX = (((Lf / Nf) * POZ) + Xf)
    DespEsfY = (((Mf / Nf) * POZ) + Yf)

    """ Calculate the radius of curvature of the reference sphere """
    SphRefRad = np.sqrt(((DespEsfX - Xf)**2.0) + ((DespEsfY - Yf)**2.0) + (POZ**2.0))

    ArrSup = SYS.SDT
    SS = []
    for ii in range(len(ArrSup)):
        ArrSup[ii].EraseVTK()
        ss = copy.deepcopy(ArrSup[ii])
        SS.append(ss)
    SS[-1].Thickness = POZ
    """ Generate the surface that is the reference sphere """
    ss = copy.deepcopy(ArrSup[ii])
    SS.append(ss)
    SS[-1].Diameter = SphRefRad
    SS[-1].Rc = -np.sign(POZ) * np.abs(SphRefRad)
    SS[-1].Thickness = 0
    SS[-1].TiltX = Tx
    SS[-1].TiltY = Ty
    SS[-1].DespX = DespEsfX
    SS[-1].DespY = DespEsfY
    SS[-1].AxisMove = 0
    SS[-1].Order = 0

    config_1 = SYS.SETUP
    SYSTEM = Kos.system(SS, config_1, build = 0)
    ###############################################################################

    POPD = []

    RAYS = Kos.raykeeper(SYSTEM)
    SYSTEM.IgnoreVignetting()
    dCos = [L[0], M[0], N[0]]

    for i in range(len(NPx)):
        pSource = [NPx[i], NPy[i], NPz[i]]
        SYSTEM.Trace(pSource, dCos, w)
        POPD.append(np.sum(SYSTEM.OP) - 2 * SYSTEM.OP[-1])
        RAYS.push()

    SYSTEM.Trace(CR_pSource, dCos, w)

    CR_OP = np.sum(SYSTEM.OP) - 2 * SYSTEM.OP[-1]
    [CR_vx, CR_vy, CR_vz] = SYSTEM.OST_XYZ[-1]

    ERX = Pupil.Cordx
    ERY = Pupil.Cordy

    ERX = np.asarray(ERX)
    ERY = np.asarray(ERY)
    Wi = np.asarray(POPD)
    Wi = -(Wi - CR_OP)
    Wi = ((Wi * 1000.0) / w)

    P2V = (np.max(Wi) - np.min(Wi))

    return P2V, Wi, ERX, ERY, SYSTEM, RAYS
