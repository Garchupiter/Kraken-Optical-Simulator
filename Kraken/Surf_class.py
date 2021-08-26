# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 16:53:41 2021

@author: JOELHERRERAVAZQUEZ
"""
import numpy as np
import pyvista as pv
import os
import sys

currentDirectory = os.getcwd()

sys.path.insert(1, currentDirectory + '/library')

from .math_shapes_class import even_asphere, aspheric__surf, extra__surf, conic__surf, axicon__surf, error_map__surf, \
    zernike__surf, zernike_expand
from .physics_class import *


class surf:
    def __init__(self):
        pass

        self.Rc = 0.
        self.Thickness = 0.00000000001
        self.Diameter = 1.0
        self.InDiameter = 0.0
        self.k = 0.0
        self.ZNK = np.zeros(36)
        self.Glass = "AIR"
        self.DespX = 0.0
        self.DespY = 0.0
        self.DespZ = 0.0
        self.TiltX = 0.0
        self.TiltY = 0.0
        self.TiltZ = 0.0
        self.Order = 0
        self.AxisMove = 1
        self.Diff_Ord = 0.0
        self.Grating_D = 0.0
        self.Grating_Angle = 0.0
        self.ShiftX = 0.0
        self.ShiftY = 0.0
        self.Mask_Type = 0  # 0 non masked, 1 apperture, 2 Obstruction
        # noinspection PyTypeChecker
        Objeto_3D = pv.Disc(center=[0.0, 0.0, 0.0], inner=0, outer=0.001, normal=(0, 0, 1), r_res=3, c_res=3)
        Mask = pv.MultiBlock()
        Mask.append(Objeto_3D)
        self.Mask_Shape = Mask
        self.Solid_3d_stl = "None"
        self.Cylinder_Rxy_Ratio = 1  # 0 for Rc_y plane 0.5 for Rc_y = 2*Rc
        self.Axicon = 0
        self.AspherData = np.zeros(200)
        self.ExtraData = np.zeros(200)
        self.Thin_Lens = 0.0
        self.Name = ""
        self.Nm_Poss = (0, 0)
        self.Note = "None"
        self.Drawing = 1
        self.Color = [0, 0, 0]
        self.Error_map = []
        self.Res = 1

        """  self.Surface_type
        0 for conic,axicon,cilinder,extra_shape, zernikes and error_map
        1 for Difraction grating
        2 for paraxial ideal lenses
        3 for STL solids
        """
        self.Surface_type = 0
        self.SURF_FUNC = []
        self.SPECIAL_SURF_FUNC = []
        self.SURF_FUNC.append(conic__surf(0, 0, 1))
        self.General_Status = self.update()

    #############################################################################
    def warning(self):
        print(" Warning, two types of surface are set at the same time,")
        print("correct this, dont set any combination que between diffraction order,")
        print(" Thin_Lens and  Solid_3d_stl.")
        print("Default value is defined to: Surface_type=0, ")
        print("    ")
        print("If it is not used, no value should be assigned")
        return

    def build_surface_function(self):
        self.SURF_FUNC = []
        self.SPECIAL_SURF_FUNC = []
        self.Surface_type = 0
        if self.Diff_Ord != 0:
            self.Surface_type = 1
            print("Surface type 1")
        if self.Thin_Lens != 0:
            self.Surface_type = 2
            print("Surface type 2")
        if self.Solid_3d_stl != "None":
            self.Surface_type = 3
            print("Surface type 3")

        if self.Diff_Ord != 0 and self.Thin_Lens != 0:
            self.warning()
            self.Surface_type = 0
        if self.Solid_3d_stl != "None" and self.Thin_Lens != 0:
            self.warning()
            self.Surface_type = 0
        if self.Solid_3d_stl != "None" and self.Diff_Ord != 0:
            self.warning()
            self.Surface_type = 0

        # --------------------------------------------------------

        if self.Surface_type == 0:  # 0 for conic,axicon,cilinder,asphere ,extra_shape, zernikes and error_map
            self.PHYSICS = snell_refraction_vector_physics()

            if np.any(self.ZNK != 0):
                NC = len(self.ZNK)
                self.Zern_pol, self.z_pow = zernike_expand(NC)

                FUNC_0 = zernike__surf(self.ZNK, self.Zern_pol, self.z_pow, self.Diameter)
                self.SURF_FUNC.append(FUNC_0)

            if np.any(self.AspherData != 0):
                FUNC_1 = aspheric__surf(self.AspherData)
                self.SURF_FUNC.append(FUNC_1)


            if self.Rc != 0:
                FUNC_2 = conic__surf(self.Rc, self.k, self.Cylinder_Rxy_Ratio)
                self.SURF_FUNC.append(FUNC_2)

            if self.Axicon != 0:
                FUNC_3 = axicon__surf(self.Cylinder_Rxy_Ratio, self.Axicon)
                self.SURF_FUNC.append(FUNC_3)

            if np.any(self.ExtraData != 0):
                FUNC_4 = extra__surf(self.ExtraData)
                self.SURF_FUNC.append(FUNC_4)

            if len(self.Error_map) != 0:
                [X, Y, Z, SPACE] = self.Error_map
                FUNC_5 = error_map__surf(X, Y, Z, SPACE)
                self.SURF_FUNC.append(FUNC_5)



            # --------------------------------------------------------
        if self.Surface_type == 1:  # 1 for Difraction grating """
            self.PHYSICS = diffraction_grating_physics()

            FUNC = conic__surf(0, self.k, self.Cylinder_Rxy_Ratio)
            self.SURF_FUNC.append(FUNC)

            # --------------------------------------------------------
        if self.Surface_type == 2:  # 2 for paraxial ideal lenses """
            self.PHYSICS = paraxial_exact_physics()

            FUNC = conic__surf(0, self.k, self.Cylinder_Rxy_Ratio)
            self.SURF_FUNC.append(FUNC)
            self.Glass = "AIR"  # Glass is forced to be AIR
            # --------------------------------------------------------
        if self.Surface_type == 3:  # 3 for STL solids """
            self.PHYSICS = snell_refraction_vector_physics()
            FUNC = conic__surf(0, self.k, self.Cylinder_Rxy_Ratio)
            self.SURF_FUNC.append(FUNC)
        # --------------------------------------------------------

    def sigma_z(self, x, y, case):

        x = x + self.ShiftX
        y = y + self.ShiftY
        Z = np.zeros_like(x)

        N_FUNC = len(self.SURF_FUNC)
        for i in range(0, N_FUNC - case):
            Z = self.SURF_FUNC[i].calculate(x, y) + Z

        return Z

    def update(self):

        self.General_Status = np.asarray([id(self.Rc),

                                          id(self.Diameter),  # 0
                                          id(self.InDiameter),  # 1
                                          id(self.k),  # 2
                                          id(self.ZNK),  # 3
                                          id(self.AspherData),  # 4
                                          id(self.ExtraData),  # 5
                                          id(self.Cylinder_Rxy_Ratio),  # 6
                                          id(self.ShiftX),  # 7
                                          id(self.ShiftY),  # 8
                                          id(self.Axicon),  # 9
                                          id(self.Thin_Lens),  # 10
                                          id(self.Error_map),  # 11

                                          id(self.Thickness),  # 12
                                          id(self.DespX),  # 13
                                          id(self.DespY),  # 14
                                          id(self.DespZ),  # 15
                                          id(self.TiltX),  # 16
                                          id(self.TiltY),  # 17
                                          id(self.TiltZ),  # 18
                                          id(self.Order),  # 19
                                          id(self.AxisMove),  # 20

                                          id(self.Glass),  # 21

                                          id(self.Diff_Ord),  # 22
                                          id(self.Grating_D),  # 23
                                          id(self.Grating_Angle),  # 24
                                          id(self.Mask_Type),  # 25
                                          id(self.Mask_Shape),  # 26
                                          id(self.Solid_3d_stl),  # 27

                                          id(self.Name),  # 28
                                          id(self.Nm_Poss),  # 29
                                          id(self.Note),  # 30
                                          id(self.Drawing),  # 31
                                          id(self.Color),  # 32
                                          id(self.Surface_type),  # 33
                                          id(self.Res)])  # 34

        return self.General_Status



    def SaveSetup(self):
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


    def RestoreSetup(self):
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


##############################################################