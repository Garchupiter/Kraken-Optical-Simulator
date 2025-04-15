## File corrected as a result of Jianhua Jiang's comments about the consistency of the results in the article

import numpy as np
import KrakenOS as Kos
import scipy
import copy
from scipy.optimize import curve_fit

























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

###############################################################################

def Phase2(Pupil):
    "Precise determination of the exit pupil "
    
    SYS = Pupil.SYSTEM
    w = Pupil.W
    RespFieldY = Pupil.FieldY
    RespFieldX = Pupil.FieldX 
    
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
        SYS.Trace([XX[0], YY[0], ZZ[0]], [LL[0], MM[0], NN[0]], w)
        [Jx, Jy, Jz] = SYS.XYZ[-1]
        [Jl, Jm, Jn] = SYS.LMN[-1]
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
    
    POZ = BestExitPupilPos(XA, YA, ZA, LA, MA, NA)
    ##############################################################################
    """ Reverse rays are traced to define the incoming rays """
    """ OJO, estos campos deben de definirse correctamente """

    Pupil.FieldX = RespFieldX    
    Pupil.FieldY = RespFieldY

    Pupil.Pattern()
    
    X, Y, Z, L, M, N = Pupil.Pattern2FieldPlus()
    
    Pupil.Ptype = "rtheta"
    
    Pupil.rad = 1
    Pupil.theta = 0
    Pupil.Pattern()
    X0, Y0, Z0, L0, M0, N0 = Pupil.Pattern2FieldPlus()
    
    Pupil.theta = 90
    Pupil.Pattern()
    X1, Y1, Z1, L1, M1, N1 = Pupil.Pattern2FieldPlus()
    
    Pupil.theta = 180
    Pupil.Pattern()
    X2, Y2, Z2, L2, M2, N2 = Pupil.Pattern2FieldPlus()
    
    Pupil.theta = 270
    Pupil.Pattern()
    X3, Y3, Z3, L3, M3, N3 = Pupil.Pattern2FieldPlus()
    
    Rays = Kos.raykeeper(SYS)
    SYS.IgnoreVignetting()
    i = 0
    
    # Chief ray
    CR_pSource = [X[i], Y[i], Z[i]]
    dCos = [L[i], M[i], N[i]]
    SYS.Trace(CR_pSource, dCos, w)
    Rays.push()
    
    pSource = [X0[i], Y0[i], Z0[i]]
    dCos = [L0[i], M0[i], N0[i]]
    SYS.Trace(pSource, dCos, w)
    Rays.push()
    
    pSource = [X1[i], Y1[i], Z1[i]]
    dCos = [L1[i], M1[i], N1[i]]
    SYS.Trace(pSource, dCos, w)
    Rays.push()
    
    pSource = [X2[i], Y2[i], Z2[i]]
    dCos = [L2[i], M2[i], N2[i]]
    SYS.Trace(pSource, dCos, w)
    Rays.push()
    
    pSource = [X3[i], Y3[i], Z3[i]]
    dCos = [L3[i], M3[i], N3[i]]
    SYS.Trace(pSource, dCos, w)
    Rays.push()
    
    # Kos.display2d(SYS, Rays, 0)
    
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
    
    NPx, NPy, NPz = find_intersections(pSource, dCos, Px, Py, Pz)
    
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
    Tx, Ty = calculate_angles_tx_ty([Lf, Mf, Nf])
    
    """ Determine the lateral displacement of the pupil with the intersection
    of the chief ray at POZ """
    DespEsfX = (((Lf / Nf) * POZ) + Xf)
    DespEsfY = (((Mf / Nf) * POZ) + Yf)
    
    """ Calculate the radius of curvature of the reference sphere """
    SphRefRad = np.sqrt(((DespEsfX - Xf)**2.0) + ((DespEsfY - Yf)**2.0) + (POZ**2.0))
    
    ArrSup = SYS.SDT
    SS = []
    for ii in range(0, len(ArrSup)):
        ArrSup[ii].EraseVTK()
        ss = copy.deepcopy(ArrSup[ii])
        SS.append(ss)
    SS[-1].Thickness = POZ
    
    """ Generate the surface that is the reference sphere """
    ss = copy.deepcopy(ArrSup[ii])
    SS.append(ss)
    SS[-1].Diameter = SphRefRad
    SS[-1].Rc = -np.sign(SS[-2].Thickness) * np.abs(SphRefRad)
      
    
    print(Tx, Ty, "----------------------------------------------")
    if Tx > 90:
        Tx = Tx - 180
    
    if Tx < -90:
        Tx = Tx + 180
    
    if Ty > 90:
        Ty = Ty - 180
    
    if Ty < -90:
        Ty = Ty + 180
    
    
    
    print(Tx, Ty, "----------------------------------------------")
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
    
    
    

    # Kos.display2d(SYSTEM, RAYS, 0)
    return (ERX, ERY, Wi, np.abs(P2V))  


def Posc(X0, Y0, L, M, N, DeltaZ):
    X1 = ((L / N) * DeltaZ) + X0
    Y1 = ((M / N) * DeltaZ) + Y0
    return X1, Y1



def Phase(PUPIL):
#############################################################
    """ Crea una pupila con las propiedades pasadas como parámetro"""


    # alpha = abcd[0]
    # betha = abcd[1]
    # gamma = abcd[2]
    # ro = abcd[3]


    SYSTEM = PUPIL.SYSTEM
    sup = PUPIL.Surf
    W = PUPIL.W
    ApType = PUPIL.ApertureType
    ApVal = PUPIL.ApertureValue
    configuracion_1 = SYSTEM.SETUP

    W = W

    """ Se calcula la pupila """
    Pup = Kos.PupilCalc(SYSTEM, sup, W, ApType, ApVal)

    [PIX, PIY, PIZ] = PUPIL.PosPupInp
    [POX, POY, POZ] = PUPIL.PosPupOutFoc

    [PIL, PIM, PIN] = PUPIL.DirPupSal
    [POL, POM, PON] = PUPIL.DirPupSal

    Pup.Samp = PUPIL.Samp
    Pup.Ptype = PUPIL.Ptype
    Pup.FieldY = PUPIL.FieldY
    Pup.FieldX = PUPIL.FieldX
    Pup.FieldType = PUPIL.FieldType

    """se generan los rayos para trazar por la pupila"""
    x, y, z, L, M, N = Pup.Pattern2Field()

    Pup.Ptype = 'chief'
    xc, yc, zc, Lc, Mc, Nc = Pup.Pattern2Field()

    pSource_0_c = [xc[0], yc[0], zc[0]]
    dCos_c = [Lc[0], Mc[0], Nc[0]]

    """ Traza el rayo principal"""
    SYSTEM.IgnoreVignetting()
    SYSTEM.Trace(pSource_0_c, dCos_c, W)

    """ Toma las coordenadas y los cosenos directores del rayo
    principal a la salida del sistema"""

    [Xcc, Ycc, Zcc] = SYSTEM.XYZ[(- 1)]
    [Lcc, Mcc, Ncc] = SYSTEM.LMN[(- 1)]


    resf_out =  np.sqrt(((POX - Xcc)**2.0) + ((POY- Ycc)**2.0) + (POZ**2.0))

# - - - - - - - - - - - - - - - - - - - - - - - - -
    ArrSup = PUPIL.SYSTEM.SDT

    SS = []
    for ii in range(0,len(ArrSup)):
        ArrSup[ii].EraseVTK()
        ss = copy.deepcopy(ArrSup[ii])
        SS.append(ss)

    DisPup = SS[0].Thickness

# - - - - - - - - - - - - - - - - - - - - - -
    """ No mover """
    SS[0].Thickness = 0.0
    SS[-1].Thickness = np.sign(POZ)*np.abs(resf_out) #+ gamma

# - - - - - - - - - - - - - - - - - - - - - -

    """ No mover """
    INP_P = Kos.surf()
    INP_P.Rc = 0
    INP_P.Thickness = DisPup
    INP_P.Diameter = PUPIL.RadPupInp * 3.0
    INP_P.Glass = 'AIR'
    INP_P.DespX = xc[0]
    INP_P.DespY = yc[0]
    INP_P.TiltX = np.rad2deg(np.arcsin((- Mc[0])))
    INP_P.TiltY = np.rad2deg(np.arcsin((Lc[0] / np.cos(np.arcsin((- Mc[0]))))))
    INP_P.AxisMove = 0
    INP_P.Order = 1

# - - - - - - - - - - - - - - - - - - - - - -

    OUT_P = Kos.surf()
    OUT_P.Rc = -np.sign(POZ)*np.abs(resf_out) #+ ro
    OUT_P.Thickness = 0
    OUT_P.Diameter = PUPIL.RadPupOut * 3
    OUT_P.Glass = 'AIR'
    OUT_P.DespX = Xcc #+ alpha
    OUT_P.DespY = Ycc #+ betha
    OUT_P.AxisMove = 0
    OUT_P.Order = 1

    Test = Kos.surf()
    Test.Diameter = PUPIL.RadPupOut * 3

# - - - - - - - - - - - - - - - - - - - - - - - - -

    SS.insert(1, INP_P)
    SS.append(OUT_P)

    SYSTEM = Kos.system(SS, configuracion_1, build = 0)

# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Se configuran los contenedores de rayos"""
    RR = Kos.raykeeper(SYSTEM)

# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Configura al sistema para ignorar los viñeteos """
    SYSTEM.IgnoreVignetting()

# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Traza los rayos por systema"""

    for i in range(0, len(x)):
        """ Traza los rayos al diametro de la pupila"""
        pSource_0 = [x[i], y[i], z[i]]
        dCos = [L[i], M[i], N[i]]
        SYSTEM.Trace(pSource_0, dCos, W)
        RR.push()

# - - - - - - - - - - - - - - - - - - - - - - - - -

    All_Rays = np.asarray(RR.OP)

    axis = 1
    AR = (np.sum(All_Rays,axis) - All_Rays[:,0]) - All_Rays[:,-1]*2


# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Traza el rayo principal"""

    SYSTEM.Trace(pSource_0_c, dCos_c, W)

    Chief = np.asarray(SYSTEM.OP)

    CH = (np.sum(Chief) - Chief[0]) - Chief[-1]*2

# - - - - - - - - - - - - - - - - - - - - - - - - -

    Wi = CH - AR

# - - - - - - - - - - - - - - - - - - - - - - - - -

    XPUP = x - xc[0]
    YPUP = y - yc[0]


    x0, y0, z0, k0, l0, mi = RR.pick(1)
    xi, yi, zi, ki, li, mi = RR.pick(-1)


    Wi = ((Wi * 1000.0) / W)




    P2V = (np.max(Wi) - np.min(Wi))

    RR.clean()
    RR.push()

    return ((YPUP / Pup.RadPupInp), (XPUP / Pup.RadPupInp), Wi, np.abs(P2V))






# ##############################################################################
# ##############################################################################
# ##############################################################################

# def PhasePlus(PUPIL, VT):
# #############################################################
#     """ Crea una pupila con las propiedades pasadas como parámetro"""


#     # alpha = abcd[0]
#     # betha = abcd[1]
#     # gamma = abcd[2]
#     # ro = abcd[3]


#     SYSTEM = PUPIL.SYSTEM
#     sup = PUPIL.Surf
#     W = PUPIL.W
#     ApType = PUPIL.ApertureType
#     ApVal = PUPIL.ApertureValue
#     configuracion_1 = SYSTEM.SETUP
#     W = W

#     """ Se calcula la pupila """
#     Pup = Kos.PupilCalc(SYSTEM, sup, W, ApType, ApVal)

#     # [PIX, PIY, PIZ] = PUPIL.PosPupInp
#     # [POX, POY, POZ] = PUPIL.PosPupOutFoc

#     [PIL, PIM, PIN] = PUPIL.DirPupSal
#     [POL, POM, PON] = PUPIL.DirPupSal

#     # Pup.Samp = PUPIL.Samp
#     # Pup.Ptype = PUPIL.Ptype
#     # Pup.FieldY = PUPIL.FieldY
#     # Pup.FieldX = PUPIL.FieldX
#     # Pup.FieldType = PUPIL.FieldType

#     """se utilizan los rayos para trazar por la pupila"""

#     [x, y, z, L, M, N, Chief_pSource, Chief_dCos, PupInp, PupOut] = VT


#     [xc, yc, zc] = Chief_pSource
#     [Lc, Mc, Nc] = Chief_dCos


#     # pSource_0_c = [xc, yc, zc]
#     # dCos_c = [Lc, Mc, Nc]


#     [PIX, PIY, PIZ] = PupInp
#     [POX, POY, POZ] = PupOut

#     """ Traza el rayo principal"""
#     SYSTEM.IgnoreVignetting()
#     SYSTEM.Trace(Chief_pSource, Chief_dCos, W)


#     """ Toma las coordenadas y los cosenos directores del rayo
#     principal a la salida del sistema"""

#     [Xcc, Ycc, Zcc] = SYSTEM.XYZ[(- 1)]
#     [Lcc, Mcc, Ncc] = SYSTEM.LMN[(- 1)]

#     # sx, sy = Posc(Xcc, Ycc, Lcc, Mcc, Ncc, POZ)
#     # print(sx, sy, "-.-.-.-.-.-.-.-.")
#     # POX = POX + sx
#     # POY = POY + sy

#     # print(POX, POY)


#     resf_out =  np.sqrt(((POX - Xcc)**2.0) + ((POY- Ycc)**2.0) + (POZ**2.0))

# # - - - - - - - - - - - - - - - - - - - - - - - - -
#     ArrSup = PUPIL.SYSTEM.SDT

#     SS = []
#     for ii in range(0,len(ArrSup)):
#         ArrSup[ii].EraseVTK()
#         ss = copy.deepcopy(ArrSup[ii])
#         SS.append(ss)

#     DisPup = SS[0].Thickness

# # - - - - - - - - - - - - - - - - - - - - - -
#     """ No mover """
#     SS[0].Thickness = 0.0
#     SS[-1].Thickness = np.sign(POZ)*np.abs(resf_out) #+ gamma

# # - - - - - - - - - - - - - - - - - - - - - -

#     """ No mover """
#     INP_P = Kos.surf()
#     INP_P.Rc = 0
#     INP_P.Thickness = DisPup
#     INP_P.Diameter = PUPIL.RadPupInp * 3.0
#     INP_P.Glass = 'AIR'
#     INP_P.DespX = xc
#     INP_P.DespY = yc
#     INP_P.TiltX = np.rad2deg(np.arcsin((- Mc)))
#     INP_P.TiltY = np.rad2deg(np.arcsin((Lc / np.cos(np.arcsin((- Mc))))))
#     INP_P.AxisMove = 0
#     INP_P.Order = 1

# # - - - - - - - - - - - - - - - - - - - - - -

#     OUT_P = Kos.surf()
#     OUT_P.Rc = -np.sign(POZ)*np.abs(resf_out)
#     OUT_P.Thickness = 0
#     OUT_P.Diameter = PUPIL.RadPupOut * 3
#     OUT_P.Glass = 'AIR'
#     OUT_P.DespX = Xcc #+ alpha
#     OUT_P.DespY = Ycc #+ betha
#     OUT_P.AxisMove = 0
#     OUT_P.Order = 1

#     Test = Kos.surf()
#     Test.Diameter = PUPIL.RadPupOut * 3

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     SS.insert(1, INP_P)
#     SS.append(OUT_P)

#     SYSTEM = Kos.system(SS, configuracion_1, build = 0)

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Se configuran los contenedores de rayos"""
#     RR = Kos.raykeeper(SYSTEM)

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Configura al sistema para ignorar los viñeteos """
#     SYSTEM.IgnoreVignetting()

# # # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Traza los rayos por systema"""

#     for i in range(0, len(x)):
#         """ Traza los rayos al diametro de la pupila"""
#         pSource_0 = [x[i], y[i], z[i]]
#         dCos = [L[i], M[i], N[i]]
#         SYSTEM.Trace(pSource_0, dCos, W)
#         RR.push()

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     All_Rays = np.asarray(RR.OP)

#     axis = 1
#     AR = (np.sum(All_Rays,axis) - All_Rays[:,0]) - All_Rays[:,-1]*2


# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Traza el rayo principal"""

#     SYSTEM.Trace(Chief_pSource, Chief_dCos, W)

#     Chief = np.asarray(SYSTEM.OP)

#     CH = (np.sum(Chief) - Chief[0]) - Chief[-1]*2

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     Wi = CH - AR

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     XPUP = x - xc
#     YPUP = y - yc


#     x0, y0, z0, k0, l0, mi = RR.pick(1)
#     xi, yi, zi, ki, li, mi = RR.pick(-1)


#     Wi = ((Wi * 1000.0) / W)




#     P2V = (np.max(Wi) - np.min(Wi))

#     RR.clean()
#     RR.push()

#     return ((YPUP / Pup.RadPupInp), (XPUP / Pup.RadPupInp), Wi, np.abs(P2V))

#     # return SYSTEM

