import os

import numpy as np
import Kraken as kn

P_Obj = kn.surf()

P_Obj.Thickness = 2.000000000000000E+003
P_Obj.Glass = "AIR"
P_Obj.Diameter = 6.796727741707513E+002 * 2.0
P_Obj.Drawing = 0

M1 = kn.surf()
M1.Rc = -6.06044E+003
M1.Thickness = -1.774190000000000E+003 + 1.853722901194000E+000
M1.k = -1.637E+000
M1.Glass = "MIRROR"
M1.Diameter = 6.63448E+002 * 2.0
M1.InDiameter = 228.6 * 2.0
M1.DespY = 0.0
M1.TiltX = 0.0000
M1.AxisMove = 1

M2 = kn.surf()
M2.Rc = -6.06044E+003
M2.Thickness = -M1.Thickness

M2.k = -3.5782E+001
M2.Glass = "MIRROR"
M2.Diameter = 2.995730651164167E+002 * 2.0
ED0 = np.zeros(20)
ED0[2] = 4.458178314555000E-018
M2.AspherData = ED0

Vertex = kn.surf()
Vertex.Thickness = 30.0
Vertex.Glass = "AIR"
Vertex.Diameter = 600.0
Vertex.Drawing = 0

Corrector_c1 = kn.surf()
Corrector_c1.Thickness = 6.6E+000
Corrector_c1.Rc = 0
Corrector_c1.Glass = "SILICASCHOTT"  # "BK7"#"LITHOSIL-Q"
Corrector_c1.Diameter = 118.0 * 2
ED1 = np.zeros(20)
ED1[0] = 2.059727552003000E-005
ED1[1] = -1.135080732384000E-009
Corrector_c1.AspherData = ED1

Corrector_c1.ExtraData = ED1

# Corrector_c1a=kn.surf()
# Corrector_c1a.Thickness=1
# Corrector_c1a.Glass="NULL"#"BK7"#"LITHOSIL-Q"

# Corrector_c1b=kn.surf()
# Corrector_c1b.Thickness=1
# Corrector_c1b.Glass="NULL"#"BK7"#"LITHOSIL-Q"


Corrector_c2 = kn.surf()
Corrector_c2.Thickness = 341.65484183207997 - 100.0
Corrector_c2.Glass = "AIR"
Corrector_c2.Diameter = 118.0 * 2.0

currentDirectory = os.getcwd()
direc = r"Prisma.stl"

objeto = kn.surf()
objeto.Diameter = 118.0 * 2.0
objeto.Solid_3d_stl = direc
objeto.Thickness = 600
objeto.Glass = "BK7"
objeto.TiltX = 55
objeto.TiltY = 0
objeto.TiltZ = 45
objeto.DespX = 0
objeto.DespY = 0
objeto.AxisMove = 0

P_Ima = kn.surf()
P_Ima.Rc = -2000
P_Ima.Thickness = 100.0
P_Ima.Glass = "BK7"
P_Ima.Diameter = 2000.0
P_Ima.Drawing = 1

P_Ima.TiltX = 45
P_Ima.AxisMove = 1

A = [P_Obj, M1, M2, Vertex, Corrector_c1, Corrector_c2, objeto, P_Ima]
# [0     1  2   3       4            5            6                7            8       9]


configuracion_1 = kn.Kraken_setup()
Telescope = kn.system(A, configuracion_1)
Rays = kn.raykeeper(Telescope)

W = 0.633
# Telescope.parax(W)

tam = 5
rad = 6.56727741707513E+002
tsis = len(A) + 2

for gg in range(0, 10):
    for j in range(-tam, tam + 1):
        # j=0
        for i in range(-tam, tam + 1):
            x_0 = (i / tam) * rad
            y_0 = (j / tam) * rad
            r = np.sqrt((x_0 * x_0) + (y_0 * y_0))
            if r < rad:
                tet = 0.0
                pSource_0 = [x_0, y_0, 0.0]
                # print("-...............")
                dCos = [0.0, np.sin(np.deg2rad(tet)), np.cos(np.deg2rad(tet))]

                W = 0.633
                Telescope.NsTrace(pSource_0, dCos, W)
                Rays.push()

kn.display3d(Telescope, Rays, 0)

print(Telescope.EFFL)
