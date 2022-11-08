
import numpy as np

class surface_tools():
    """surface_tools.
    """
    def __init__(self, SurfData):
        """__init__.

        Parameters
        ----------
        SurfData :
            SurfData
        """
        self.SDT = SurfData
        self.ErrSurfCase = 0
        self.Surface_Flattener = (- 1)

    def SurfaceShape(self, x, y, j):
        """SurfaceShape.

        Parameters
        ----------
        x :
            x
        y :
            y
        j :
            j
        """
        if (j != self.Surface_Flattener):
            TOTAL_SURF_SHAPE = self.SDT[j].sigma_z(x, y, self.ErrSurfCase)
        else:
            TOTAL_SURF_SHAPE = 0.0*np.copy(x)
        return TOTAL_SURF_SHAPE

