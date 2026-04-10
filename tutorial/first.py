"""Examp Doublet Lens Pupil"""

# Loading the library
import KrakenOS as Kos

# Creating an object of the surf class for the object plane
P_Obj = Kos.surf()
P_Obj.Thickness = 100
P_Obj.Glass = "AIR"
P_Obj.Diameter = 30.0

# Creating a surface for the first face in BK7 Glass
L1a = Kos.surf()
L1a.Rc = 92.847
L1a.Thickness = 6.0
L1a.Glass = "BK7"
L1a.Diameter = 30.0

# Creating a surface for the second face in F2 Glass
L1b = Kos.surf()
L1b.Rc = -30.716
L1b.Thickness = 3.0
L1b.Glass = "F2"
L1b.Diameter = 30

# Creating a surface for the third interface to air
L1c = Kos.surf()
L1c.Rc = -78.197
L1c.Thickness = 97.376 - 40
L1c.Glass = "AIR"
L1c.Diameter = 30

# Creating a surface to exemplify a pupil
pupil = Kos.surf()
pupil.Rc = 30
pupil.Thickness = 40.
pupil.Glass = "AIR"
pupil.Diameter = 5
pupil.Name = "Pupil"
pupil.DespY = 0.
pupil.Nm_Poss=[-10,10]

# Creating a surface for image plane
P_Ima = Kos.surf()
P_Ima.Rc = 0.0
P_Ima.Thickness = 0.0
P_Ima.Glass = "AIR"
P_Ima.Diameter = 20.0
P_Ima.Name = "P_Ima"
P_Ima.Nm_Poss=[-10,10]

A = [P_Obj, L1a, L1b, L1c, pupil, P_Ima]
config_1 = Kos.Setup()

# Creating the system with previus information
Doublet = Kos.system(A, config_1)

Rays = Kos.raykeeper(Doublet)

W = 0.4
sur = 4
AperVal = 10
AperType = "EPD"
Pup = Kos.PupilCalc(Doublet, sur, W, AperType, AperVal)

# Configuring field and ray array type
Pup.Samp = 3
Pup.Ptype = "fan"
Pup.FieldType = "angle"
Pup.FieldY = 2.0


# ray origin coordinates and direction cosines
x, y, z, L, M, N = Pup.Pattern2Field()

# Tracing the rays with a loop
for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Doublet.Trace(pSource_0, dCos, W)
    Rays.push()# Saving rays

# Configuring (-field) and ray array type,.. etc
Pup.FieldY = -Pup.FieldY
x, y, z, L, M, N = Pup.Pattern2Field()
for i in range(0, len(x)):
    pSource_0 = [x[i], y[i], z[i]]
    dCos = [L[i], M[i], N[i]]
    Doublet.Trace(pSource_0, dCos, W)
    Rays.push() # Saving rays

Kos.display2d(Doublet, Rays,0,1)
