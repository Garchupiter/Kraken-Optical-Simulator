
import numpy as np
import pyvista as pv

class Prerequisites():
    """Prerequisites.
    """

    def __init__(self, SurfData, SUTO):
        """__init__.

        Parameters
        ----------
        SurfData :
            SurfData
        SUTO :
            SUTO
        """
        self.SDT = SurfData
        self.n = len(self.SDT)
        self.ANG = np.arange(0, 362, 1)
        self.ANG_Rad = np.deg2rad(self.ANG)
        self.Disable_Inner = 1
        self.ExtraDiameter = 0
        self.SuTo = SUTO
        self.TRANS_1A = []
        self.TRANS_2A = []

    def GeometricRotatAndTran(self, L_te_h, j):
        """GeometricRotatAndTran.

        Parameters
        ----------
        L_te_h :
            L_te_h
        j :
            j
        """
        stx = 0
        sty = 0
        for n in range(j, (- 1), (- 1)):
            LL = 1.0
            if (n == 0):
                LL = 0
            if (n == j):
                PA = 1
            else:
                PA = (self.SDT[n].AxisMove * LL)
            tx = ((self.SDT[n].TiltX * PA) * LL)
            ty = ((self.SDT[n].TiltY * PA) * LL)
            tz = ((self.SDT[n].TiltZ * PA) * LL)
            dx = ((self.SDT[n].DespX * PA) * LL)
            dy = ((self.SDT[n].DespY * PA) * LL)
            dz = ((self.SDT[n].DespZ * PA) * LL)
            stx = (stx + tx)
            sty = (sty + ty)
            Tol_Err = 0.0001
            if (np.abs(np.cos(np.deg2rad(stx))) < Tol_Err):
                tx = (tx + Tol_Err)
            if (np.abs(np.cos(np.deg2rad(sty))) < Tol_Err):
                ty = (ty + Tol_Err)
            if (self.SDT[n].Order == 0):
                L_te_h.rotate_x(tx, inplace=True)
                L_te_h.rotate_y(ty, inplace=True)
                L_te_h.rotate_z(-tz, inplace=True)
                L_te_h.translate([dx, dy, dz],inplace=True)
                L_te_h.translate([0, 0, self.SDT[(n - 1)].Thickness],inplace=True)
            else:
                L_te_h.translate([dx, dy, dz],inplace=True)
                L_te_h.rotate_z(-tz,inplace=True)
                L_te_h.rotate_y(ty,inplace=True)
                L_te_h.rotate_x(tx,inplace=True)
                L_te_h.translate([0, 0, self.SDT[(n - 1)].Thickness],inplace=True)
        return L_te_h

    def Flat2SigmaSurface(self, plane_object, j):
        """Flat2SigmaSurface.

        Parameters
        ----------
        plane_object :
            plane_object
        j :
            j
        """
        plane_objectAx = plane_object.points[:, 0]
        plane_objectAy = plane_object.points[:, 1]
        plane_objectAz = plane_object.points[:, 2]

        plane_objectAx = np.asarray(plane_objectAx)
        plane_objectAy = np.asarray(plane_objectAy)
        plane_objectAz = np.asarray(plane_objectAz)

        plane_objectAx = plane_objectAx + self.SDT[j].SubAperture[2]
        plane_objectAy = plane_objectAy + self.SDT[j].SubAperture[1]

        plane_objectAz = self.SuTo.SurfaceShape(plane_objectAx, plane_objectAy, j)
        plane_objectC = np.c_[(plane_objectAx, plane_objectAy, plane_objectAz)]
        plane_object.points = plane_objectC
        return plane_object

    def Face3D(self, j):
        """Face3D.

        Parameters
        ----------
        j :
            j
        """
        if self.SDT[j].Mask_Shape == "None":
            self.SDT[j].RestoreVTK()


        if (self.SDT[j].Solid_3d_stl == 'None'):
            RES = (47 * self.SDT[j].Res)
            con = (((self.SDT[j].Diameter *self.SDT[j].SubAperture[0]) - (self.SDT[j].InDiameter * self.Disable_Inner)) / (self.SDT[j].Diameter *self.SDT[j].SubAperture[0]))
            if (con == 0):
                con = 1
            r_RES = int((RES * con))
            INNER = ((self.SDT[j].InDiameter * self.Disable_Inner) / 2.0)
            OUTER = ((self.SDT[j].Diameter *self.SDT[j].SubAperture[0]) / 2.0)
            disc = pv.Disc(center=[0.0, 0.0, 0.0], inner=INNER, outer=OUTER, normal=(0, 0, 1), r_res=r_RES, c_res=(RES * 2))
            L_te_h = self.Flat2SigmaSurface(disc, j)
            if (self.SDT[j].InDiameter > 0):
                L_te_h = L_te_h.delaunay_2d().edge_source = L_te_h
            else:
                L_te_h = L_te_h.delaunay_2d()
        else:
            aaa = isinstance(self.SDT[j].Solid_3d_stl, pv.core.pointset.PolyData)
            if aaa == False:
                L_te_h = pv.read(self.SDT[j].Solid_3d_stl)
            else:
                L_te_h = self.SDT[j].Solid_3d_stl

            L_te_h.compute_normals(cell_normals=True, point_normals=True, split_vertices=True, flip_normals=False, consistent_normals=True, auto_orient_normals=False, non_manifold_traversal=True, feature_angle=30.0, inplace = False)
        L_te_h = self.GeometricRotatAndTran(L_te_h, j)
        MASK = self.SDT[j].Mask_Shape
        OBJECT_MASK = pv.MultiBlock()
        for mask in MASK:
            Mask_poly = self.Flat2SigmaSurface(mask, j)
            Mask_poly = self.GeometricRotatAndTran(Mask_poly, j)
            OBJECT_MASK.append(Mask_poly)

        return (L_te_h, OBJECT_MASK)

    def SidePerim(self, j):
        """SidePerim.

        Parameters
        ----------
        j :
            j
        """
        rad2 = ((self.SDT[j].Diameter *self.SDT[j].SubAperture[0]) / 2.0)

        x2 = rad2 * np.cos(self.ANG_Rad)
        y2 = rad2 * np.sin(self.ANG_Rad)
        z2 = np.zeros_like(x2)
        x2 = x2 + self.SDT[j].SubAperture[2]
        y2 = y2 + self.SDT[j].SubAperture[1]

        points2 = np.c_[(x2, y2, z2)]
        L_te = pv.PolyData(points2, force_float=False)
        x2 = L_te.points[:, 0]
        y2 = L_te.points[:, 1]
        z2 = self.SuTo.SurfaceShape(x2, y2, j)
        points2 = np.c_[(x2, y2, z2)]
        L_te = pv.PolyData(points2, force_float=False)
        L_te = self.GeometricRotatAndTran(L_te, j)
        return L_te.points

    def Side3D(self, j, j2):
        """Side3D.

        Parameters
        ----------
        j :
            j
        j2 :
            j2
        """
        PTS1 = self.SidePerim(j)
        x1 = PTS1[:, 0]
        y1 = PTS1[:, 1]
        z1 = PTS1[:, 2]
        PTS2 = self.SidePerim(j2)
        x2 = PTS2[:, 0]
        y2 = PTS2[:, 1]
        z2 = PTS2[:, 2]

        RR = np.zeros_like(x2)
        for s in range(0, len(x2)):
            x2 = np.roll(x2, 1)
            y2 = np.roll(y2, 1)
            z2 = np.roll(z2, 1)

            X = x2 - x1
            Y = y2 - y1
            Z = z2 - z1
            R = np.sqrt(X**2 + Y**2 + Z**2)


            RR[s] = np.mean(R)
        ag = np.argmin(RR)

        x2 = np.roll(x2, ag)
        y2 = np.roll(y2, ag)
        z2 = np.roll(z2, ag)

        P = []
        F = []
        for i in range(0, 360):
            P.append([x1[i], y1[i], z1[i]])
            P.append([x2[i], y2[i], z2[i]])
            P.append([x1[(i + 1)], y1[(i + 1)], z1[(i + 1)]])
            F.append([3, (0 + (i * 6)), (1 + (i * 6)), (2 + (i * 6))])
            P.append([x1[(i + 1)], y1[(i + 1)], z1[(i + 1)]])
            P.append([x2[i], y2[i], z2[i]])
            P.append([x2[(i + 1)], y2[(i + 1)], z2[(i + 1)]])
            F.append([3, (3 + (i * 6)), (4 + (i * 6)), (5 + (i * 6))])
        P = np.asarray(P)
        F = np.asarray(F)
        cant = pv.PolyData(P, F, force_float=False)
        Ax = cant.points[:, 0]
        Ay = cant.points[:, 1]
        Az = cant.points[:, 2]
        Ax = np.asarray(Ax)
        Ay = np.asarray(Ay)
        Az = np.asarray(Az)

        cant.compute_normals(cell_normals=True, point_normals=True, split_vertices=True, flip_normals=False, consistent_normals=True, auto_orient_normals=False, non_manifold_traversal=True, feature_angle=30.0, inplace = False)

        return cant

    def Prerequisites3DSolids(self):
        """Prerequisites3DSolids.
        """
        self.GlassOnSide = []
        self.PreTypeTotal = []
        self.TypeTotal = []
        self.PreGlassOnSide = []
        self.AAA = pv.MultiBlock()
        self.BBB = pv.MultiBlock()
        self.DDD = pv.MultiBlock()
        self.EEE = pv.MultiBlock()
        self.counter = []
        self.side_number = []
        for j in range(0, self.n):
            (lens, masked) = self.Face3D(j)
            self.AAA.append(lens)
            if (self.SDT[j].Solid_3d_stl == 'None'):
                self.TypeTotal.append(0)
            else:
                self.TypeTotal.append(1)
            self.EEE.append(lens)
            self.counter.append(j)
            self.DDD.append(masked.copy(deep=True))
            self.counter.append(j)
            self.GlassOnSide.append(j)
            if (j < (self.n - 1)):
                if (self.SDT[j].Glass != 'NULL'):
                    if (self.SDT[j].Glass != 'AIR'):
                        if (self.SDT[j].Glass != 'MIRROR'):
                            if (self.SDT[j].Solid_3d_stl == 'None'):
                                j2 = (j + 1)
                                while True:
                                    if ((self.SDT[j].Glass == 'NULL') or (self.SDT[j2].Solid_3d_stl != 'None')):
                                        j2 = (j2 + 1)
                                    else:
                                        break
                                    if (j2 == (self.n - 1)):
                                        break
                                side = self.Side3D(j, j2)
                                self.BBB.append(side)
                                self.PreTypeTotal.append(1)
                                self.PreGlassOnSide.append(j)
                                self.side_number.append(j)
        for i in self.BBB:
            self.EEE.append(i)
        for i in self.PreGlassOnSide:
            self.GlassOnSide.append(i)
        for i in self.PreTypeTotal:
            self.TypeTotal.append(i)

    def Prerequisites3SMath(self):
        """Prerequisites3SMath.
        """
        self.TRANS_1A = []
        self.TRANS_2A = []
        for j in range(0, self.n):
            Start_trans1 = np.matrix([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]])
            for n in range(j, (- 1), (- 1)):
                LL = 1.0
                if (n == 0):
                    LL = 0.0
                if (n == j):
                    PA = 1
                else:
                    PA = (self.SDT[n].AxisMove * LL)
                Tx = (((- np.deg2rad(self.SDT[n].TiltX)) * PA) * LL)
                Ty = (((- np.deg2rad(self.SDT[n].TiltY)) * PA) * LL)
                Tz = (((- np.deg2rad(self.SDT[n].TiltZ)) * PA) * LL)
                TH = self.SDT[(n - 1)].Thickness
                (dx, dy, dz) = (((self.SDT[n].DespX * PA) * LL), ((self.SDT[n].DespY * PA) * LL), ((self.SDT[n].DespZ * PA) * LL))
                DTH_Z = np.matrix([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, (- TH)], [0.0, 0.0, 0.0, 1.0]])
                Rx_1A = np.matrix([[1.0, 0.0, 0.0, 0.0], [0.0, np.cos((- Tx)), np.sin((- Tx)), 0.0], [0.0, (- np.sin((- Tx))), np.cos((- Tx)), 0.0], [0.0, 0.0, 0.0, 1.0]])
                Ry_1A = np.matrix([[np.cos((- Ty)), 0.0, (- np.sin((- Ty))), 0.0], [0, 1.0, 0.0, 0.0], [np.sin((- Ty)), 0.0, np.cos((- Ty)), 0.0], [0.0, 0.0, 0.0, 1.0]])
                Rz_1A = np.matrix([[np.cos((- Tz)), (- np.sin((- Tz))), 0.0, 0.0], [np.sin((- Tz)), np.cos((- Tz)), 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0, 0, 0, 1.0]])
                Dxyz_1A = np.matrix([[1.0, 0.0, 0.0, (- dx)], [0.0, 1.0, 0.0, (- dy)], [0.0, 0.0, 1.0, (- dz)], [0.0, 0.0, 0.0, 1.0]])
                if (self.SDT[n].Order == 0):
                    Start_trans1 = (((((Start_trans1 * Rz_1A) * Ry_1A) * Rx_1A) * Dxyz_1A) * DTH_Z)
                else:
                    Start_trans1 = (((((Dxyz_1A * Rx_1A) * Ry_1A) * Rz_1A) * Start_trans1) * DTH_Z)
            self.TRANS_1A.append(Start_trans1)
            self.TRANS_2A.append(np.linalg.inv(Start_trans1))

