import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import scipy
import pkg_resources
required = {'KrakenOS'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    print("No instalado")
    import sys
    sys.path.append("../..")


import KrakenOS as Kos



# ______________________________________#

def CosDirsSums(v1, v2):

    [L1, M1, N1] = v1

    y1 =  (M1 / N1)
    x1 =  (L1 / N1)

    tetx1 = np.arctan2(x1, 1.0)
    tety1 = np.arctan2(y1, 1.0)

    [L2, M2, N2] = v2

    y2 =  (M2 / N2)
    x2 =  (L2 / N2)

    tetx2 = np.arctan2(x2, 1.0)
    tety2 = np.arctan2(y2, 1.0)

    tetx = tetx1 + tetx2
    tety = tety1 + tety2

    x2 = np.tan(tetx)
    y2 = np.tan(tety)
    z2 = 1.0
    r = np.sqrt((x2**2) + (y2**2) + (z2**2))

    L = x2/r
    M = y2/r
    N = z2/r

    LMN = np.asarray([L, M, N])
    return(LMN)

# ______________________________________#

P_Obj = Kos.surf()
P_Obj.Rc = 0
P_Obj.Thickness = 1000 + 3.452200000000000E+003
P_Obj.Glass = "AIR"
P_Obj.Diameter = 1.059E+003 * 2.0

# ______________________________________#

Thickness = 3.452200000000000E+003
M1 = Kos.surf()
M1.Rc = -9.638000000004009E+003
M1.Thickness = -Thickness
M1.k = -1.077310000000000E+000
M1.Glass = "MIRROR"
M1.Diameter = 1.059E+003 * 2.0
M1.InDiameter = 250 * 2.0
M1.TiltY = 0.0
M1.TiltX = 0.0

# ______________________________________#

M1.AxisMove = 0
M2 = Kos.surf()
M2.Rc = -3.93E+003
M2.Thickness = Thickness + 1037.525880
M2.k = -4.328100000000000E+000
M2.Glass = "MIRROR"
M2.Diameter = 3.365E+002 * 2.0
M2.TiltY = 0.0
M2.TiltX = 0.0
M2.DespY = 0.0
M2.DespX = 0.0
M2.AxisMove = 0

# ______________________________________#

P_Ima = Kos.surf()
P_Ima.Diameter = 300.0
P_Ima.Glass = "AIR"
P_Ima.Name = "Plano imagen"

# ______________________________________#

A = [P_Obj, M1, M2, P_Ima]
configuracion_1 = Kos.Setup()
Telescope = Kos.system(A, configuracion_1)
Rays1 = Kos.raykeeper(Telescope)
Rays2 = Kos.raykeeper(Telescope)



Surf = 1
W = 0.5016
AperVal = 2000.
AperType = "EPD"
Pupil = Kos.PupilCalc(Telescope, Surf, W, AperType, AperVal)
Pupil.Samp = 20
Pupil.Ptype = "rand"
Pupil.FieldX = 0.00025
Pupil.FieldY = 0.0
Pupil.FieldType = "angle"

# Estrella 1
x1, y1, z1, L1, M1, N1 = Pupil.Pattern2Field()

Pupil.FieldX = -0.00025
# Estrella 2
x2, y2, z2, L2, M2, N2 = Pupil.Pattern2Field()


# ------------------------------------------------#


Sun = Kos.SourceRnd()

# Gaussian (Seeing)
def f(x):
    x = np.rad2deg(x)
    seing = 1.2 / 3600.0
    sigma = seing / 2.3548
    mean = 0
    standard_deviation = sigma
    y = scipy.stats.norm(mean, standard_deviation)
    res = y.pdf(x)
    return res


Sun.field = 4 * 1.2 / (2.0 * 3600.0)

Sun.fun = f
Sun.dim = 2300

# ------------------------------------------------------

Sun.num = len(x1)
L0, M0, N0, X0, Y0, Z0 = Sun.rays()

v1 = np.array([L0, M0, N0])

v2 = np.array([L1, M1, N1])

[L, M, N] = CosDirsSums(v1, v2)


for i in range(0, len(x1)):
    pSource_0 = [x1[i], y1[i], z1[i]]
    dCos = [L[i], M[i], N[i]]
    Telescope.Trace(pSource_0, dCos, W)
    Rays1.push()

# ------------------------------------------------------


Sun.num = len(x2)
L0, M0, N0, X0, Y0, Z0 = Sun.rays()

v1 = np.array([L0, M0, N0])

v2 = np.array([L2, M2, N2])

[L, M, N] = CosDirsSums(v1, v2)


for i in range(0, len(x2)):
    pSource_0 = [x2[i], y2[i], z2[i]]
    dCos = [L[i], M[i], N[i]]
    Telescope.Trace(pSource_0, dCos, W)
    Rays2.push()
# ------------------------------------------------------



# Kos.display3d(Telescope, Rays, 1)


X1, Y1, Z1, L1, M1, N1 = Rays1.pick(-1)
X2, Y2, Z2, L2, M2, N2 = Rays2.pick(-1)

# ______________________________________#




plt.plot(X1, Y1, 'x', c = "b")
plt.plot(X2, Y2, 'o', c = "r")

plt.xlabel('numbers')
plt.ylabel('values')
plt.title('spot Diagram')
plt.axis('square')
plt.show()


print("Hola mundo")
