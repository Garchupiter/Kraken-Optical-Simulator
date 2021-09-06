# Kraken-Optical-Simulator (KrakenOS)          
## Python - Exact ray tracing library 

Joel Herrera V., Carlos Guerrero P., Morgan Rhaí Najera Roa, Anais Sotelo B., Ilse Plauchu F.

• joel@astro.unam.mx

KrakenOS (Kraken - Optical Simulator) is a python library based in Numpy, Matplotlib, PyVTK and PyVista libraries, it provides a three-dimensional optical systems visualization and ray tracing. This tool has been programed on the object-oriented paradigm. KrakenOS focuses on performing sequential and non-sequential exact ray tracing, it permits to define all the parameters of the optical elements or even the mathematical function to describe their shape, it also allows adding optical properties to 3D solid elements in STL format and use glass catalogs. The library permit to control and modifying the position of the surfaces in a three-dimensional space, this allows generating off-axis systems. It also has several tools such as the calculation of wavefront aberrations in terms of Zernike polynomials, Seidel sums, Entrance and exit pupil calculation and paraxial optics.

## Prerequisites
The library has been tested with the following packages and versions.
• Python '3.7.4'          
• numpy '1.18.5'          
• pyvista '0.25.3'          
• pyvtk '0.5.18'  
• matplotlib '3.4.3'  
• vtk '8.2'          
• csv '1.0'          
• Place the directory “KrakenOS” in the same path where the code to be executed is located.          

## Surfaces and the optical system
The library has been simplified to the point of having only two classes of objects for the definition of a system, these are surf and system.
The surf object contains all the relevant information of every optical interface, in this way, every optical interface is an object of the surf class, all interfaces, from the object plane to the image plane, contain attributes of size, shape, material or orien##tation.

## User manual and examples
Very important, read the included user manual (KrakenOS_User_Manual.pdf) and the set of useful examples:
• Examp_Axicon.py          
• Examp_Axicon_And_Cylinder.py          
• Examp_Diffraction_Grating_Reflection.py          
• Examp_Diffraction_Grating_Transmission.py          
• Examp_Doublet_Lens_ParaxMatrix.py          
• Examp_Doublet_Lens.py          
• Examp_Doublet_Lens_3Dcolor.py          
• Examp_Doublet_Lens_CommandsSystem.py          
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
