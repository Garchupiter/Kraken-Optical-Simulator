# Kraken-Optical-Simulator
Python - Exact ray tracing library

Joel Herrera V., Carlos Guerrero P., Morgan Rhaí Najera Roa, Anais Sotelo B., Ilse Plauchu F.          

• joel@astro.unam.mx           


Kraken consists of a library developed in the Python 3.0 using Numpy, MAtplotlib, PyVTK and PyVista libraries which allows it to provide a three-dimensional visualization of the optical elements. This tool has focused on the paradigm of object-oriented programming , this seems to be naturally the one that gives us greater simplicity in the implementation of a system. Kraken focuses on performing the exact ray tracing in a sequential and non-sequential way, it also allows to vary all the parameters of the optical elements or even define the mathematical function to describe their shape, it also allows adding optical properties to 3D solid elements in format STL and use glass catalogs. Another of the library's capabilities is that it allows modifying the position of surfaces in a three-dimensional space, this allows generating off-axis systems, it has several tools such as the calculation of wavefront aberrations in terms of Zernike polynomials, Seidel sums, Entrance and exit pupil calculation and paraxial optics.
## Prerequisites

The library has been tested with the following packages and versions.
• Python '3.7.4'          
• numpy '1.18.5'          
• pyvista '0.25.3'          
• pyvtk '0.5.18'  
• matplotlib '3.4.3'  
• vtk '8.2'          
• csv '1.0'          
• Place the directory “Kraken” in the same path where the code to be executed is located.          

## Surfaces and the optical system
The library has been simplified to the point of having only two classes of objects for the definition of a system, these are surf and system.
The surf object contains all the relevant information of every optical interface, in this way, every optical interface is an object of the surf class, all interfaces, from the object plane to the image plane, contain attributes of size, shape, material or orien##tation.

## User manual and examples
Very important, read the user manual (Kraken_User_Manual.pdf) included and the set of useful examples included with the library.

• Examp-Axicon.py          
• Examp-Axicon_And_Cylinder.py          
• Examp-Diffraction_Grating_Reflection.py          
• Examp-Diffraction_Grating_Transmission.py          
• Examp-Doublet_Lens-ParaxMatrix.py          
• Examp-Doublet_Lens.py          
• Examp-Doublet_Lens_3Dcolor.py          
• Examp-Doublet_Lens_CommandsSystem.py          
• Examp-Doublet_Lens_Cylinder.py          
• Examp-Doublet_Lens_NonSec.py          
• Examp-Doublet_Lens_Pupil.py          
• Examp-Doublet_Lens_Pupil_Seidel.py          
• Examp-Doublet_Lens_Tilt-Nulls.py          
• Examp-Doublet_Lens_Tilt.py          
• Examp-Doublet_Lens_Tilt_non_sec.py          
• Examp-Doublet_Lens_Zernike.py          
• Examp-ExtraShape_Micro_Lens_Array.py          
• Examp-ExtraShape_Radial_Sine.py          
• Examp-ExtraShape_XY_Cosines.py          
• Examp-Flat_Mirror_45Deg.py          
• Examp-MultiCore.py          
• Examp-ParaboleMirror_Shift.py          
• Examp-Perfect_lens.py          
• Examp-Ray.py          
• Examp-Solid_Object_STL.py          
• Examp-Solid_Object_STL_ARRAY.py          
• Examp-Sourse_Distribution_Function.py          
• Examp-Tel_2M.py          
• Examp-Tel_2M_Error_Map.py          
• Examp-Tel_2M_Pupila.py          
• Examp-Tel_2M_Spyder_Spot_Diagram.py          
• Examp-Tel_2M_Spyder_Spot_RMS.py          
• Examp-Tel_2M_Spyder_Spot_Tilt_M2.py          
• Examp-Tel_2M_Atmospheric_Refraction.py          
• Examp-Tel_2M_Wavefront_Fitting.py         

Enjoy it!!
