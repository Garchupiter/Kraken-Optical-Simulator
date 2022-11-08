
import numpy as np
import pyvista as pv
import os
import sys
currentDirectory = os.getcwd()
sys.path.insert(1, (currentDirectory + '/library'))
from .MathShapesClass import *
from .PhysicsClass import *

class surf():
    """surf.
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
             Ratio between the axial and sagittal
             radius of curvature.


          surf.Axicon = 0
             values other than zero an axicon is
             generated with the angle defined


          surf.Thickness = 0.0
             Distance between this surface and the next surface.


          surf.Diameter = 1.0
             Outside diameter of the surface.


          surf.InDiameter = 0.0
             Internal diameter of the surface.


          surf.k = 0.0
             Conicity constant for classical conic surfaces
             k = 0 for spherical
             k = -1 for parabola, etc.
             Default value: 0.0


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
             Defines what will happen to the optical axis
             after a coordinate transformation.


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
             (0) Do not apply mask,
             (1) Use mask as aperture,
             (2) Use mask as obstruction.
             Default value: 0


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

        """


    def __init__(self, Rc = 0.0, Thickness = 1e-11,Diameter = 1.0, InDiameter = 0.0, k = 0.0, Glass = 'AIR'\
        ,ZNK = np.zeros(36)\
        ,DespX = 0.0\
        ,DespY = 0.0\
        ,DespZ = 0.0\
        ,TiltX = 0.0\
        ,TiltY = 0.0\
        ,TiltZ = 0.0\
        ,Order = 0.0\
        ,AxisMove = 1.0\
        ,Diff_Ord = 0.0\
        ,Grating_D = 0.0\
        ,Grating_Angle = 0.0\
        ,ShiftX = 0.0\
        ,ShiftY = 0.0\
        ,Mask_Type = 0.0\

        ,Mask_Shape = "None"\
        ,Solid_3d_stl = "None"\
        ,Cylinder_Rxy_Ratio = 1.0\
        ,Axicon = 0.0\
        ,AspherData = np.zeros(200)\
        ,ExtraData = np.zeros(200)\
        ,Thin_Lens = 0.0\
        ,Name = ''\
        ,Nm_Pos = (0.0, 0.0)\
        ,Note = "None"\
        ,Drawing = 1.0\
        ,Color = [0, 0, 0]\
        ,Error_map = []\
        ,Res = 1\
        ,Surface_type = 0.0\
        ,SURF_FUNC = [conic__surf(0.0, 0.0, 1.0)]\
        ,SPECIAL_SURF_FUNC = []\
        ,SubAperture = [1,0,0], Coating = [[],[],[],[]], NumLabel = 1):

        """__init__.
        """
        pass
        self.Rc = Rc
        self.Thickness = Thickness
        self.Diameter = Diameter
        self.InDiameter = InDiameter
        self.k = k
        self.Glass = Glass

        self.ZNK = ZNK
        self.DespX = DespX
        self.DespY = DespY
        self.DespZ = DespZ
        self.TiltX = TiltX
        self.TiltY = TiltY
        self.TiltZ = TiltZ
        self.Order = Order
        self.AxisMove = AxisMove
        self.Diff_Ord = Diff_Ord
        self.Grating_D = Grating_D
        self.Grating_Angle = Grating_Angle
        self.ShiftX = ShiftX
        self.ShiftY = ShiftY
        self.Mask_Type = Mask_Type

        self.Mask_Shape = Mask_Shape
        self.Solid_3d_stl = Solid_3d_stl
        self.Cylinder_Rxy_Ratio = Cylinder_Rxy_Ratio
        self.Axicon = Axicon
        self.AspherData = AspherData
        self.ExtraData = ExtraData
        self.Thin_Lens = Thin_Lens
        self.Name = Name
        self.Nm_Pos = Nm_Pos
        self.Note = Note
        self.Drawing = Drawing
        self.Color = Color
        self.Error_map = Error_map
        self.Res = Res
        self.Surface_type = Surface_type
        self.SURF_FUNC = SURF_FUNC
        self.SPECIAL_SURF_FUNC = SPECIAL_SURF_FUNC
        self.General_Status = self.update()
        self.SubAperture = SubAperture
        self.Coating = Coating
        self.NumLabel = NumLabel

    def RestoreVTK(self):
        Objeto_3D = pv.Disc(center=[0.0, 0.0, 0.0], inner=0, outer=0.001, normal=(0, 0, 1), r_res=3, c_res=3)
        Mask = pv.MultiBlock()
        Mask.append(Objeto_3D)
        self.Mask_Shape = Mask

    def EraseVTK(self):

        self.Mask_Shape = "None"

    def warning(self):
        """warning.
        """
        print(' Warning, two types of surface are set at the same time,')
        print('correct this, dont set any combination que between diffraction order,')
        print(' Thin_Lens and  Solid_3d_stl.')
        print('Default value is defined to: Surface_type=0, ')
        print('    ')
        print('If it is not used, no value should be assigned')
        return

    def build_surface_function(self):
        """build_surface_function.
        """
        self.SURF_FUNC = []
        self.SPECIAL_SURF_FUNC = []
        self.Surface_type = 0
        if (self.Diff_Ord != 0):
            self.Surface_type = 1
            # print('Surface type 1')
        if (self.Thin_Lens != 0):
            self.Surface_type = 2
            # print('Surface type 2')
        if (self.Solid_3d_stl != 'None'):
            self.Surface_type = 3
            # print('Surface type 3')
        if ((self.Diff_Ord != 0) and (self.Thin_Lens != 0)):
            self.warning()
            self.Surface_type = 0
        if ((self.Solid_3d_stl != 'None') and (self.Thin_Lens != 0)):
            self.warning()
            self.Surface_type = 0
        if ((self.Solid_3d_stl != 'None') and (self.Diff_Ord != 0)):
            self.warning()
            self.Surface_type = 0
        if (self.Surface_type == 0):
            self.PHYSICS = snell_refraction_vector_physics()
            if np.any((self.ZNK != 0)):
                NC = len(self.ZNK)
                (self.Zern_pol, self.z_pow) = zernike_expand(NC)
                FUNC_0 = zernike__surf(self.ZNK, self.Zern_pol, self.z_pow, self.Diameter)
                self.SURF_FUNC.append(FUNC_0)
            if np.any((self.AspherData != 0)):
                FUNC_1 = aspheric__surf(self.AspherData)
                self.SURF_FUNC.append(FUNC_1)
            if (self.Rc != 0):
                FUNC_2 = conic__surf(self.Rc, self.k, self.Cylinder_Rxy_Ratio)
                self.SURF_FUNC.append(FUNC_2)
            if (self.Axicon != 0):
                FUNC_3 = axicon__surf(self.Cylinder_Rxy_Ratio, self.Axicon)
                self.SURF_FUNC.append(FUNC_3)
            if np.any((self.ExtraData != 0)):
                FUNC_4 = extra__surf(self.ExtraData)
                self.SURF_FUNC.append(FUNC_4)
            if (len(self.Error_map) != 0):
                [X, Y, Z, SPACE] = self.Error_map
                FUNC_5 = error_map__surf(X, Y, Z, SPACE)
                self.SURF_FUNC.append(FUNC_5)
        if (self.Surface_type == 1):
            self.PHYSICS = diffraction_grating_physics()
            FUNC = conic__surf(0, self.k, self.Cylinder_Rxy_Ratio)
            self.SURF_FUNC.append(FUNC)
        if (self.Surface_type == 2):
            self.PHYSICS = paraxial_exact_physics()
            FUNC = conic__surf(0, self.k, self.Cylinder_Rxy_Ratio)
            self.SURF_FUNC.append(FUNC)
            self.Glass = 'AIR'
        if (self.Surface_type == 3):
            self.PHYSICS = snell_refraction_vector_physics()
            FUNC = conic__surf(0, self.k, self.Cylinder_Rxy_Ratio)
            self.SURF_FUNC.append(FUNC)

    def sigma_z(self, x, y, case):
        """sigma_z.

        Parameters
        ----------
        x :
            x
        y :
            y
        case :
            case
        """
        x = (x + self.ShiftX)
        y = (y + self.ShiftY)

        con = -1
        N_FUNC = len(self.SURF_FUNC)

        for i in range(0, (N_FUNC - case)):

            if con==-1:
                Z = self.SURF_FUNC[i].calculate(x, y)
                con = 0

            else:
                Z = (self.SURF_FUNC[i].calculate(x, y) + Z)
                con=1

        if con == -1:
            Z = 0.0 * np.copy(x)

        return Z



    def update(self):
        """update.
        """
        self.General_Status = np.asarray([id(self.Rc), id(self.Diameter), id(self.InDiameter), id(self.k), id(self.ZNK), id(self.AspherData), id(self.ExtraData), id(self.Cylinder_Rxy_Ratio), id(self.ShiftX), id(self.ShiftY), id(self.Axicon), id(self.Thin_Lens), id(self.Error_map), id(self.Thickness), id(self.DespX), id(self.DespY), id(self.DespZ), id(self.TiltX), id(self.TiltY), id(self.TiltZ), id(self.Order), id(self.AxisMove), id(self.Glass), id(self.Diff_Ord), id(self.Grating_D), id(self.Grating_Angle), id(self.Mask_Type), id(self.Mask_Shape), id(self.Solid_3d_stl), id(self.Name), id(self.Nm_Pos), id(self.Note), id(self.Drawing), id(self.Color), id(self.Surface_type), id(self.Res)])
        return self.General_Status

    def SaveSetup(self):
        """SaveSetup.
        """
        self.Sv_Rc = self.Rc
        self.Sv_Thickness = self.Thickness
        self.Sv_Diameter = self.Diameter
        self.Sv_InDiameter = self.InDiameter
        self.Sv_k = self.k
        self.Sv_ZNK = self.ZNK
        self.Sv_Glass = self.Glass
        self.Sv_DespX = self.DespX
        self.Sv_DespY = self.DespY
        self.Sv_DespZ = self.DespZ
        self.Sv_TiltX = self.TiltX
        self.Sv_TiltY = self.TiltY
        self.Sv_TiltZ = self.TiltZ
        self.Sv_Order = self.Order
        self.Sv_AxisMove = self.AxisMove
        self.Sv_Diff_Ord = self.Diff_Ord
        self.Sv_Grating_D = self.Grating_D
        self.Sv_Grating_Angle = self.Grating_Angle
        self.Sv_ShiftX = self.ShiftX
        self.Sv_ShiftY = self.ShiftY
        self.Sv_Mask_Type = self.Mask_Type
        self.Sv_Cylinder_Rxy_Ratio = self.Cylinder_Rxy_Ratio
        self.Sv_Axicon = self.Axicon
        self.Sv_AspherData = self.AspherData
        self.Sv_ExtraData = self.ExtraData
        self.Sv_Thin_Lens = self.Thin_Lens
        self.Sv_Error_map = self.Error_map
        self.Sv_SubAperture = self.SubAperture
        self.Sv_Coating = self.Coating

    def RestoreSetup(self):
        """RestoreSetup.
        """
        self.Rc = self.Sv_Rc
        self.Thickness = self.Sv_Thickness
        self.Diameter = self.Sv_Diameter
        self.InDiameter = self.Sv_InDiameter
        self.k = self.Sv_k
        self.ZNK = self.Sv_ZNK
        self.Glass = self.Sv_Glass
        self.DespX = self.Sv_DespX
        self.DespY = self.Sv_DespY
        self.DespZ = self.Sv_DespZ
        self.TiltX = self.Sv_TiltX
        self.TiltY = self.Sv_TiltY
        self.TiltZ = self.Sv_TiltZ
        self.Order = self.Sv_Order
        self.AxisMove = self.Sv_AxisMove
        self.Diff_Ord = self.Sv_Diff_Ord
        self.Grating_D = self.Sv_Grating_D
        self.Grating_Angle = self.Sv_Grating_Angle
        self.ShiftX = self.Sv_ShiftX
        self.ShiftY = self.Sv_ShiftY
        self.Mask_Type = self.Sv_Mask_Type
        self.Cylinder_Rxy_Ratio = self.Sv_Cylinder_Rxy_Ratio
        self.Axicon = self.Sv_Axicon
        self.AspherData = self.Sv_AspherData
        self.ExtraData = self.Sv_ExtraData
        self.Thin_Lens = self.Sv_Thin_Lens
        self.Error_map = self.Sv_Error_map
        self.SubAperture = self.Sv_SubAperture
        self.Coating = self.Sv_Coating
