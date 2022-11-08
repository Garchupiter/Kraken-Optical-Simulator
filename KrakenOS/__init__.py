#!/usr/bin/env Python3
"""
 Kraken-Optical-Simulator (KrakenOS)
 Python - Exact ray tracing library

Joel Herrera V., Carlos Guerrero P., Morgan Rhaí Najera Roa, Anais Sotelo B., Ilse Plauchu F.

• joel@astro.unam.mx

KrakenOS (Kraken - Optical Simulator) is a python library based in Numpy, Matplotlib, PyVTK and PyVista libraries, it provides a three-dimensional optical systems visualization and ray tracing. This tool has been programed on the object-oriented paradigm. KrakenOS focuses on performing sequential and non-sequential exact ray tracing, it permits to define all the parameters of the optical elements or even the mathematical function to describe their shape, it also allows adding optical properties to 3D solid elements in STL format and use glass catalogs. The library permit to control and modifying the position of the surfaces in a three-dimensional space, this allows generating off-axis systems. It also has several tools such as the calculation of wavefront aberrations in terms of Zernike polynomials, Seidel sums, Entrance and exit pupil calculation and paraxial optics.

PREREQUISITES
The library has been tested with the following packages and versions.
• Python '3.7.4'
• numpy '1.18.5'
• pyvista '0.25.3'
• pyvtk '0.5.18'
• matplotlib '3.4.3'
• vtk '8.2'
• csv '1.0'
• Place the directory “KrakenOS” in the same path where the code to be executed is located.

Surfaces and the optical system
The library has been simplified to the point of having only two classes of objects for the definition of a system, these are surf and •   system.
The surf object contains all the relevant information of every optical interface, in this way, every optical interface is an object of the surf class, all interfaces, from the object plane to the image plane, contain attributes of size, shape, material or orientation.

SURF CLASS ATRIBUTES:

  surf.Name = ""
     Name of the element.

  surf.NamePos = (0,0)
     “Name” position in the 2D diagram.

  surf.Note = "None"
     Useful for adding user notes to a surface.

  surf.Rc = 0
     Paraxial radius of curvature in millimeters.

  surf.Cylinder_Rxy_Ratio = 1
     Ratio between the axial and sagittal radius of curvature.

  surf.Axicon = 0
     values other than zero an axicon is generated with the angle defined

  surf.Thickness = 0.0
     Distance between this surface and the next surface.

  surf.Diameter = 1.0
     Outside diameter of the surface.

  surf.InDiameter = 0.0
     Internal diameter of the surface.

  surf.k = 0.0
     Conicity constant for classical conic surfaces, k = 0 for spherical, k = -1 for parabola, etc. Default value: 0.0

  surf.DespX = 0.0
  surf.DespY = 0.0
  surf.DespZ = 0.0
    Displacement of the surface in the X, Y and Z axis

  surf.TiltX = 0.0
  surf.TiltY = 0.0
  surf.TiltZ = 0.0
     Rotation of the surface in the X, Y and Z axis

  surf.Order = 0
     Define the order of the transformations.

  surf.AxisMove = 1
     Defines what will happen to the optical axis after a coordinate transformation.

  surf.Diff_Ord = 0.0
     Diffraction order.

  surf.Grating_D = 0.0
     Separation between the lines of the diffraction grating.

  surf.Grating_Angle = 0.0
     Angle of the grating lines in the plane of the surface

  surf.ZNK = np.zeros ()
     Zernike polynomials coefficients

  surf.ShiftX = 0
  surf.ShiftY = 0
     Offset the surface function on the X or Y axis.

  surf.Mask = 0
     (0) Do not apply mask, (1) Use mask as aperture, (2) Use mask as obstruction. Default value: 0

  surf.Mask_Shape = Object_3D
     Form of the mask to apply on surface

  surf.AspherData = np.zeros ()
     Asphericity coefficients.

  surf.ExtraData = [f, coef]
     User-defined function for optical surface

  surf.Error_map = [X, Y, Z, SPACE]
     Error map array

  surf.Drawing = 1
     1 for drawn in the 3D plot, 0 to omit.

  surf.Color = [0,0,0]
     Element color for 3D Plot. [R,G,B]

  surf.Solid_3d_stl = "None"
     Path to the 3D solid STL file.


SYSTEM CLASS ATRIBUTES:

  system.Trace (pS, dC, wV)
     Sequential ray tracing.
     pS = [1.0, 0.0, 0.0] – Ray origin coordinates
     dC = [0.0,0.0,1.0] - The directing cosines
     wV = 0.4 - Wavelength

  system.NsTrace (pS, dC, wV)
     Non-Sequential ray tracing

  system.Parax (w)
     Paraxial optics calculations

  system.disable_inner
     Enables the central aperture.

  system.enable_inner
	Disables the central aperture.

  system.SURFACE
	Returns the surfaces the ray passed through.

  system.NAME
	Returns surface names that the ray passed through

  system.GLASS
	Returns materials that the ray passed through.

  system.XYZ
	[X, Y, Z] ray coordinates from its origin to the image plane.

  system.OST_XYZ
	[X, Y, Z] coordinates of ray intersections with respect to a coordinate system at its vertex, even if this vertex has a translation or rotation.

  system.DISTANCE
     List of distances traveled by the ray.

  system.OP
	List of optical paths.

  system.TOP
	Total optical path.

  system.TOP_S     List of the ray's optical path by sections.

  system.ALPHA
     List the materials absorption coefficients

  system.BULK_TRANS
     List the transmission through all the •   system. absorption coefficients are considered.

  system.S_LMN
     Surface normal direction cosines [L, M, N].

  system.LMN
     incident ray direction cosines [L, M, N].

  system.R_LMN
     Resulting ray direction cosines [L, M, N].

  system.N0
     Refractive indices before and after each interface

  system.N1
     Refractive indices after each interface. This is useful to differentiate between index before and after an iteration. Example:
     N0 = [n1, n2, n3, n4, n5]  N1 = [n2, n3, n4, n5, n5]

  system.WAV
	Wavelength of the ray (µm)

  system.G_LMN
 	[L, M, N] Direction cosines that define the lines of the diffraction grating on the plane.

  system.ORDER
	      Ray diffraction order.

  system.GRATING_D
	 Distance between lines of the diffraction grating.Units (µm)

  system.RP
	Fresnel reflection coefficients for polarization P.

  system.RS
	Fresnel reflection coefficients for polarization S.

  system.TP
	Fresnel transmission coefficients for polarization P.

  system.TS
	Fresnel transmission coefficients for polarization S.

  system.TTBE
	Total energy transmitted or reflected per element.

  system.TT
	Total energy transmitted or reflected total.

  system.targ_surf (int)
    Limits the ray tracing to the defined surface

  system.flat_surf (int)
    Change a surface to flat.



USER MANUAL AND EXAMPLES
Very important, read the included user manual (KrakenOS_User_Manual.pdf) and the set of useful examples:

• Examp_Axicon.py
• Examp_Axicon_And_Cylinder.py
• Examp_Diffraction_Grating_Reflection.py
• Examp_Diffraction_Grating_Transmission.py
• Examp_Doublet_Lens_ParaxMatrix.py
• Examp_Doublet_Lens.py
• Examp_Doublet_Lens_3Dcolor.py
• Examp_Doublet_Lens_Commands•   system.py
• Examp_Doublet_Lens_Cylinder.py
• Examp_Doublet_Lens_NonSec.py
• Examp_Doublet_Lens_Pupil.py
• Examp_Doublet_Lens_Pupil_Seidel.py
• Examp_Doublet_Lens_Tilt_Nulls.py
• Examp_Doublet_Lens_Tilt.py
• Examp_Doublet_Lens_Tilt_non_sec.py
• Examp_Doublet_Lens_Zernike.py
• Examp_ExtraShape_Micro_Lens_Array.py
• Examp_ExtraShape_Radial_Sine.py
• Examp_ExtraShape_XY_Cosines.py
• Examp_Flat_Mirror_45Deg.py
• Examp_MultiCore.py
• Examp_ParaboleMirror_Shift.py
• Examp_Perfect_lens.py
• Examp_Ray.py
• Examp_Solid_Object_STL.py
• Examp_Solid_Object_STL_ARRAY.py
• Examp_Sourse_Distribution_Function.py
• Examp_Tel_2M.py
• Examp_Tel_2M_Error_Map.py
• Examp_Tel_2M_Pupila.py
• Examp_Tel_2M_Spyder_Spot_Diagram.py
• Examp_Tel_2M_Spyder_Spot_RMS.py
• Examp_Tel_2M_Spyder_Spot_Tilt_M2.py
• Examp_Tel_2M_Atmospheric_Refraction.py
• Examp_Tel_2M_Wavefront_Fitting.py

Enjoy it!!

"""


# import inspect

from .Physics import *
from .Display import *
from .PhysicsClass import *
from .SurfClass import *
from .PupilTool import *
from .RayKeeper import *
from .SystemTools import *
from .SetupClass import *
from .SurfTools import surface_tools as SUT
from .Prerequisites3D import *
from .HitOnSurf import *
from .InterNormalCalc import *
from .WavefrontFit import *
from .SeidelTool import *
from .SourceRand import *
from .LibRMS import *
from .TraceLoopTool import *
from .PhaseCalc import *
from .WavePlot import *
from .KrakenSys import *
from .PSFCalc import *
