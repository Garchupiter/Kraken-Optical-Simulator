# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 17:17:14 2021

@author: JOELHERRERAVAZQUEZ
"""
import numpy as np

class surface_tools:

    def __init__(self, SurfData):
        """surface_tools."""

        """__init__.

        :param SurfData:
        """
        self.SDT = SurfData

        self.ErrSurfCase = 0
        self.Surface_Flattener = -1

    def SurfaceShape(self, x, y, j):
        """SurfaceShape.

        :param x:
        :param y:
        :param j:
        """
        if j != self.Surface_Flattener:
            self.Surface_Flattener

            TOTAL_SURF_SHAPE = self.SDT[j].sigma_z(x, y, self.ErrSurfCase)
        else:
            TOTAL_SURF_SHAPE = np.zeros_like(x)
        return TOTAL_SURF_SHAPE
