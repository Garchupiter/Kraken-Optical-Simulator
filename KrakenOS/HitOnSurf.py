# from numba import jit
# @jit(forceobj=True)

import numpy as np
from .SurfTools import surface_tools as SUT
from scipy.optimize import fsolve
import math

class Hit_Solver():
    """Hit_Solver.
    """


    def __init__(self, SurfData):
        """__init__.

        Parameters
        ----------
        SurfData :
            SurfData
        """
        self.SDT = SurfData
        self.n = len(self.SDT)
        self.SuTo = SUT(SurfData)
        self.vj = 0
        self.vevaX = 0
        self.vevaY = 0
        self.NP_x1 = 0
        self.NP_y1 = 0
        self.NP_z1 = 0
        self.MN = 0
        self.LN = 0
        self.SuTo.ErrSurfCase = 0
        self.delta = 0.05
        self.delta2 = 2.0 * self.delta
        self.deltaLineCurve = 1e-13
        self.Parray=np.asarray([0.0, 0.0, 0.0])

        self.vx = np.asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.vy = np.asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.vz = np.asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        self.vxi = [self.delta, 0.0, 0.0, - self.delta, 0.0, 0.0]
        self.vyi = [0.0, self.delta, 0.0, 0.0, - self.delta, 0.0]
        self.vzi = [0.0, 0.0, self.delta, 0.0, 0.0, - self.delta]



    # def SurfDer(self, x, y, z):
    #     """SurfDer.

    #     Parameters
    #     ----------
    #     x :
    #         x
    #     y :
    #         y
    #     z :
    #         z
    #     """

    #     # self.vx = np.asarray([(x + self.delta), x, x, (x - self.delta), x, x])
    #     # self.vy = np.asarray([y, (y + self.delta), y, y, (y - self.delta), y])
    #     # self.vz = np.asarray([z, z, (z + self.delta), z, z, (z - self.delta)])


    #     self.vx1 = np.asarray([(x + self.delta), x, x, (x - self.delta), x, x, (x + self.delta2), x, x, (x - self.delta2), x, x])
    #     self.vy1 = np.asarray([y, (y + self.delta), y, y, (y - self.delta), y, y, (y + self.delta2), y, y, (y - self.delta2), y])
    #     self.vz1 = np.asarray([z, z, (z + self.delta), z, z, (z - self.delta), z, z, (z + self.delta2), z, z, (z - self.delta2)])


    #     F1 = self.SuTo.SurfaceShape(self.vx1, self.vy1, self.vj) - self.vz1


    #     V = F1 / (2.0 * self.delta)


    #     Dx = V[0] - V[3]
    #     Dy = V[1] - V[4]
    #     Dz = V[2] - V[5]

    #     Dr = np.sqrt((((Dx ** 2.) + (Dy ** 2.)) + (Dz ** 2.)))

    #     return ((Dx / Dr), (Dy / Dr), (Dz / Dr))

    def SurfDer(self, x, y, z):
        """SurfDer.

        Parameters
        ----------
        x :
            x
        y :
            y
        z :
            z
        """



        self.vx1 = np.asarray([(x + self.delta), x, x, (x - self.delta), x, x, (x + self.delta2), x, x, (x - self.delta2), x, x])
        self.vy1 = np.asarray([y, (y + self.delta), y, y, (y - self.delta), y, y, (y + self.delta2), y, y, (y - self.delta2), y])
        self.vz1 = np.asarray([z, z, (z + self.delta), z, z, (z - self.delta), z, z, (z + self.delta2), z, z, (z - self.delta2)])


        F = self.SuTo.SurfaceShape(self.vx1, self.vy1, self.vj) - self.vz1

        AX1=-F[6] #-(self.__xyzF(x+(2*delta),y,z,j))
        AX2=8.*F[0] #8.0*(self.__xyzF(x+delta,y,z,j))
        AX3=-8.*F[3] #-8.0*(self.__xyzF(x-delta,y,z,j))
        AX4=F[9] #(self.__xyzF(x-(2*delta),y,z,j) )
        Dx=(AX1+AX2+AX3+AX4)/(12.0*self.delta)


        AY1=-F[6+1] #-(self.__xyzF(x+(2*delta),y,z,j))
        AY2=8.*F[0+1] #8.0*(self.__xyzF(x+delta,y,z,j))
        AY3=-8.*F[3+1] #-8.0*(self.__xyzF(x-delta,y,z,j))
        AY4=F[9+1] #(self.__xyzF(x-(2*delta),y,z,j) )
        Dy=(AY1+AY2+AY3+AY4)/(12.0*self.delta)



        AZ1=-F[6+2] #-(self.__xyzF(x+(2*delta),y,z,j))
        AZ2=8.*F[0+2] #8.0*(self.__xyzF(x+delta,y,z,j))
        AZ3=-8.*F[3+2] #-8.0*(self.__xyzF(x-delta,y,z,j))
        AZ4=F[9+2] #(self.__xyzF(x-(2*delta),y,z,j) )
        Dz=(AZ1+AZ2+AZ3+AZ4)/(12.0*self.delta)



        # V = F1 / (2.0 * self.delta)


        # Dx = V[0] - V[3]
        # Dy = V[1] - V[4]
        # Dz = V[2] - V[5]

        Dr = np.sqrt((((Dx ** 2.) + (Dy ** 2.)) + (Dz ** 2.)))

        return ((Dx / Dr), (Dy / Dr), (Dz / Dr))








    def __surface_Derivative(self,x,y,z,j):
        #http://www2.math.umd.edu/~dlevy/classes/amsc466/lecture-notes/differentiation-chap.pdf
        delta=0.000001
        # more precision but slower execution
        AX1=-(self.__xyzF(x+(2*delta),y,z,j))
        AX2=8.0*(self.__xyzF(x+delta,y,z,j))
        AX3=-8.0*(self.__xyzF(x-delta,y,z,j))
        AX4=(self.__xyzF(x-(2*delta),y,z,j) )
        Dx=(AX1+AX2+AX3+AX4)/(12*delta)

        AY1=-(self.__xyzF(x,y+2*delta,z,j))
        AY2=8.0*(self.__xyzF(x, y+delta,z,j))
        AY3=-8.0*(self.__xyzF(x,y-delta,z,j))
        AY4=(self.__xyzF(x, y-2*delta,z,j)  )
        Dy=(AY1+AY2+AY3+AY4)/(12*delta)

        AZ1=-(self.__xyzF(x,y,z+2*delta,j))
        AZ2=8.0*(self.__xyzF(x,y, z+delta,j))
        AZ3=-8.0*(self.__xyzF(x,y,z-delta,j))
        AZ4=(self.__xyzF(x, y, z-2*delta,j) )
        Dz=(AZ1+AZ2+AZ3+AZ4)/(12*delta)

        # Dx=(self.__xyzF(x+delta,y,z,j)-self.__xyzF(x-delta,y,z,j))/(2.0*delta)
        # Dy=(self.__xyzF(x,y+delta,z,j)-self.__xyzF(x,y-delta,z,j))/(2.0*delta)
        # Dz=(self.__xyzF(x,y,z+delta,j)-self.__xyzF(x,y,z-delta,j))/(2.0*delta)

        Dr=np.sqrt((Dx*Dx)+(Dy*Dy)+(Dz*Dz))
        return Dx/Dr,Dy/Dr,Dz/Dr








    def __DerLineCurve(self, gf):
        """__DerLineCurve.

        Parameters
        ----------
        gf :
            gf
        """
        self.Parray[0]= gf + self.deltaLineCurve
        self.Parray[1]= gf
        self.Parray[2]= gf - self.deltaLineCurve

        kl=self.Parray - self.NP_z1
        X = (( kl * self.LN ) + self.NP_x1)
        Y = (( kl * self.MN ) + self.NP_y1)


        SP = self.SuTo.SurfaceShape(X, Y, self.vj)
        AA = (SP - self.Parray)

        Dz = ((AA[0] - AA[2]) / (2.0 * self.deltaLineCurve))
        return Dz, AA[1]

    def SolveHit(self, Px1, Py1, Pz1, L, M, N, j):
        """SolveHit.

        Parameters
        ----------
        Px1 :
            Px1
        Py1 :
            Py1
        Pz1 :
            Pz1
        L :
            L
        M :
            M
        N :
            N
        j :
            j
        """

        # if (len(self.SDT[j].Error_map) == 0):
        #     self.SuTo.ErrSurfCase = 0
        # else:
        #     self.SuTo.ErrSurfCase = 0

        self.vevaX = 0.
        self.vevaY = 0.
        self.vj = j

        # cnt = 0
        PP_z2 = Pz1
        P_z2 = Pz1
        self.NP_x1 = Px1
        self.NP_y1 = Py1
        self.NP_z1 = Pz1
        self.MN = (M / N)
        self.LN = (L / N)

        for i in range (0,30):

            DerFdeX, FdeX = self.__DerLineCurve(PP_z2)

            PP_z2 = (PP_z2 - (FdeX / DerFdeX))

            if (np.abs((PP_z2 - P_z2)) > 1e-9):
                P_z2 = PP_z2
            else:
                break



        P_z2 = PP_z2
        P_x2 = (((P_z2 - self.NP_z1) * self.LN) + self.NP_x1)
        P_y2 = (((P_z2 - self.NP_z1) * self.MN) + self.NP_y1)

        # self.SuTo.ErrSurfCase = 0


        if (len(self.SDT[j].Error_map) != 0):
            PP_z2 = P_z2
            cnt = 0
            while True:
                DerFdeX , FdeX= self.__DerLineCurve(PP_z2)
                PP_z2  = (PP_z2 - (FdeX / DerFdeX))
                if (np.abs((PP_z2 - P_z2)) < 1e-9):
                    break
                else:
                    P_z2 = PP_z2
                if (cnt == 3):
                    break
                cnt += 1
            P_z2 = PP_z2
            P_x2 = (((P_z2 - self.NP_z1) * self.LN) + self.NP_x1)
            P_y2 = (((P_z2 - self.NP_z1) * self.MN) + self.NP_y1)

        self.vevaX = P_x2
        self.vevaY = P_y2


        return (P_x2, P_y2, P_z2)















# import numpy as np
# from .SurfTools import surface_tools as SUT
# from scipy.optimize import fsolve
# import math

# class Hit_Solver():
#     """Hit_Solver.
#     """


#     def __init__(self, SurfData):
#         """__init__.

#         Parameters
#         ----------
#         SurfData :
#             SurfData
#         """
#         self.SDT = SurfData
#         self.n = len(self.SDT)
#         self.SuTo = SUT(SurfData)
#         self.vj = 0
#         self.vevaX = 0
#         self.vevaY = 0
#         self.NP_x1 = 0
#         self.NP_y1 = 0
#         self.NP_z1 = 0
#         self.MN = 0
#         self.LN = 0
#         self.SuTo.ErrSurfCase = 0
#         self.delta = 0.0001
#         self.deltaLineCurve = 1e-13
#         self.Parray=np.asarray([0., 0.])


#     def __Fxyz(self, x, y, z):
#         """__Fxyz.

#         Parameters
#         ----------
#         x :
#             x
#         y :
#             y
#         z :
#             z
#         """

#         F = self.SuTo.SurfaceShape(x, y, self.vj)

#         VF = F - z
#         return VF

#     def SurfDer(self, x, y, z):
#         """SurfDer.

#         Parameters
#         ----------
#         x :
#             x
#         y :
#             y
#         z :
#             z
#         """
#         # print(np.shape(x), np.shape(y), np.shape(z),  "ññlk ñlk ñlk ñl ")

#         vx = np.asarray([(x + self.delta), x, x, (x - self.delta), x, x])
#         vy = np.asarray([y, (y + self.delta), y, y, (y - self.delta), y])
#         vz = np.asarray([z, z, (z + self.delta), z, z, (z - self.delta)])


#         V = self.__Fxyz(vx, vy, vz)/ (2.0 * self.delta)


#         Dx = V[0] - V[3]
#         Dy = V[1] - V[4]
#         Dz = V[2] - V[5]

#         Dr = np.sqrt((((Dx ** 2.) + (Dy ** 2.)) + (Dz ** 2.)))

#         return ((Dx / Dr), (Dy / Dr), (Dz / Dr))

#     def __lineCurve(self, vz):
#         """__lineCurve.

#         Parameters
#         ----------
#         vz :
#             vz
#         """
#         kl=vz - self.NP_z1
#         self.vevaX = ((kl * self.LN) + self.NP_x1)
#         self.vevaY = ((kl * self.MN) + self.NP_y1)


#         Vzf = (self.SuTo.SurfaceShape(self.vevaX, self.vevaY, self.vj) - vz)

#         return Vzf

#     def __DerLineCurve(self, gf):
#         """__DerLineCurve.

#         Parameters
#         ----------
#         gf :
#             gf
#         """
#         self.Parray[0]= gf + self.deltaLineCurve
#         self.Parray[1]= gf - self.deltaLineCurve

#         kl=self.Parray - self.NP_z1
#         X = (( kl * self.LN ) + self.NP_x1)
#         Y = (( kl * self.MN ) + self.NP_y1)

#         AA = (self.SuTo.SurfaceShape(X, Y, self.vj) - self.Parray)

#         Dz = ((AA[0] - AA[1]) / (2.0 * self.deltaLineCurve))
#         return Dz

#     def SolveHit(self, Px1, Py1, Pz1, L, M, N, j):
#         """SolveHit.

#         Parameters
#         ----------
#         Px1 :
#             Px1
#         Py1 :
#             Py1
#         Pz1 :
#             Pz1
#         L :
#             L
#         M :
#             M
#         N :
#             N
#         j :
#             j
#         """
#         # if (len(self.SDT[j].Error_map) != 0):
#         #     self.SuTo.ErrSurfCase = 1

#         if (len(self.SDT[j].Error_map) == 0):
#             self.SuTo.ErrSurfCase = 0
#         else:
#             self.SuTo.ErrSurfCase = 0

#         self.vevaX = 0.
#         self.vevaY = 0.
#         self.vj = j

#         # cnt = 0
#         PP_z2 = Pz1
#         P_z2 = Pz1
#         self.NP_x1 = Px1
#         self.NP_y1 = Py1
#         self.NP_z1 = Pz1
#         self.MN = (M / N)
#         self.LN = (L / N)

#         for i in range (0,30):
#             FdeX = self.__lineCurve(PP_z2)

#             DerFdeX = self.__DerLineCurve(PP_z2)
#             PP_z2 = (PP_z2 - (FdeX / DerFdeX))

#             if (np.abs((PP_z2 - P_z2)) > 1e-05):
#                 P_z2 = PP_z2
#             else:
#                 break


#             # if (np.abs((PP_z2 - P_z2)) < 1e-05):
#             #     break
#             # else:
#             #     P_z2 = PP_z2
#             # if (cnt == 30):
#             #     break
#             # cnt += 1


#         P_z2 = PP_z2
#         P_x2 = (((P_z2 - self.NP_z1) * self.LN) + self.NP_x1)
#         P_y2 = (((P_z2 - self.NP_z1) * self.MN) + self.NP_y1)

#         self.vevaX = P_x2
#         self.vevaY = P_y2
#         self.SuTo.ErrSurfCase = 0


#         if (len(self.SDT[j].Error_map) != 0):
#             PP_z2 = P_z2
#             cnt = 0
#             while True:
#                 FdeX = self.__lineCurve(PP_z2)
#                 DerFdeX = self.__DerLineCurve(PP_z2)
#                 PP_z2 = (PP_z2 - (FdeX / DerFdeX))
#                 if (np.abs((PP_z2 - P_z2)) < 1e-05):
#                     break
#                 else:
#                     P_z2 = PP_z2
#                 if (cnt == 3):
#                     break
#                 cnt += 1
#             P_z2 = PP_z2
#             P_x2 = (((P_z2 - self.NP_z1) * self.LN) + self.NP_x1)
#             P_y2 = (((P_z2 - self.NP_z1) * self.MN) + self.NP_y1)
#         return (P_x2, P_y2, P_z2)

