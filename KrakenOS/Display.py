
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
import sys
from matplotlib import rc


def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = 0
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    Px = xdata[end_ind] - (xdata[end_ind] - xdata[start_ind]) / 2.0
    Py = ydata[end_ind] - (ydata[end_ind] - ydata[start_ind]) / 2.0
    PPx = xdata[start_ind] + (xdata[end_ind] - xdata[start_ind]) / 2.01
    PPy = ydata[start_ind] + (ydata[end_ind] - ydata[start_ind]) / 2.01
    line.axes.annotate('',
        xytext=(PPx,PPy),
        xy=(Px, Py),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size)

###############################################################################

def wavelength_to_rgb(wavelength, gamma=1.0):
    """wavelength_to_rgb.

    Parameters
    ----------
    wavelength :
        wavelength
    gamma :
        gamma
    """
    wavelength = float(wavelength)
    if (380 <= wavelength <= 440):
        attenuation = (0.3 + ((0.7 * (wavelength - 380)) / (440 - 380)))
        R = ((((- (wavelength - 440)) / (440 - 380)) * attenuation) ** gamma)
        G = 0.0
        B = ((1.0 * attenuation) ** gamma)
    elif (440 <= wavelength <= 490):
        R = 0.0
        G = (((wavelength - 440) / (490 - 440)) ** gamma)
        B = 1.0
    elif (490 <= wavelength <= 510):
        R = 0.0
        G = 1.0
        B = (((- (wavelength - 510)) / (510 - 490)) ** gamma)
    elif (510 <= wavelength <= 580):
        R = (((wavelength - 510) / (580 - 510)) ** gamma)
        G = 1.0
        B = 0.0
    elif (580 <= wavelength <= 645):
        R = 1.0
        G = (((- (wavelength - 645)) / (645 - 580)) ** gamma)
        B = 0.0
    elif (645 <= wavelength <= 750):
        attenuation = (0.3 + ((0.7 * (750 - wavelength)) / (750 - 645)))
        R = ((1.0 * attenuation) ** gamma)
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    R *= 255
    G *= 255
    B *= 255
    return [(int(R) / 255.0), (int(G) / 255.0), (int(B) / 255.0)]

###############################################################################

def display3d_old(SYSTEM, RAYS, view=0, inline=False,     BackgCol= 'white', BackgColTop = 'white', GridCol="black"):

    """display3d.

    Parameters
    ----------
    SYSTEM :
        SYSTEM
    RAYS :
        RAYS
    view :
        view
    """

    ST1="KrakenOS v0.1. Executing Script: " + sys.argv[0]
    OPA = 0.95

    CCC = pv.MultiBlock()
    INST = isinstance(RAYS, list)
    if INST == True:
        for R in RAYS:
            for rays in R.CC:
                RAY_VTK_OBJ = pv.lines_from_points(rays)
                CCC.append(RAY_VTK_OBJ)
    else:
        for rays in RAYS.CC:
            RAY_VTK_OBJ = pv.lines_from_points(rays)
            CCC.append(RAY_VTK_OBJ)


    p = pv.Plotter(shape=(1, 1), title=ST1,notebook=inline)
    Absorb_color = np.array([(10 / 256), (23 / 256), (24 / 256)])
    Mirror_color = np.array([(189 / 256), (189 / 256), (189 / 256)])
    Glass_color=np.array([12/256, 238/256, 246/256])
    recorte = view
    NN = SYSTEM.AAA.n_blocks
    if (SYSTEM.SDT[0].Drawing == 0):
        points2 = np.c_[(0, 0, 0)]
        c = pv.PolyData(points2)
        cc = pv.PolyData(points2)
    if (SYSTEM.SDT[0].Drawing == 1):
        c = SYSTEM.AAA[0]
        cc = SYSTEM.AAA[0]
    for n in range(1, NN):
        if (SYSTEM.SDT[n].Drawing == 1):
            AAAA = SYSTEM.AAA[n]
            if (SYSTEM.SDT[n].Glass != 'NULL'):
                if (SYSTEM.SDT[n].Color == [0, 0, 0]):
                    if (SYSTEM.SDT[n].Glass == 'MIRROR'):
                        color = Mirror_color
                    else:
                        color = Glass_color
                    if (SYSTEM.SDT[n].Glass == 'ABSORB'):
                        color = Absorb_color
                else:
                    color = SYSTEM.SDT[n].Color
                cc = cc.merge(AAAA)
                if (recorte == 1):
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    clippedy = clippedx.clip((0, 1, 0), invert=False)
                    c = c.merge(clippedy)
                    clippedx = AAAA.clip(((- 1), 0, 0), invert=False)
                    clippedy = clippedx.clip((0, (- 1), 0), invert=False)
                    c = c.merge(clippedy)
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    clippedy = clippedx.clip((0, (- 1), 0), invert=False)
                    c = c.merge(clippedy)
                if (recorte == 0):
                    c = cc
                if (recorte == 2):
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    c = c.merge(clippedx)
                p.add_mesh(c, color, opacity=OPA, specular=1, specular_power=15, smooth_shading=True, show_edges=False)
                edges = c.extract_feature_edges(feature_angle=10, boundary_edges=True, feature_edges=False, manifold_edges=False)
                if SYSTEM.SDT[n].Solid_3d_stl != "None" and recorte ==0:
                    print(" ") # No edges
                else:
                    p.add_mesh(edges, 'red')
                points2 = np.c_[(0, 0, 0)]
                c = pv.PolyData(points2)
    p.add_mesh(SYSTEM.DDD, color=[0.5, 0.5, 0.5], opacity=OPA, show_edges=None)
    NN = SYSTEM.AAA.n_blocks
    n = 0
    for g in SYSTEM.side_number:
        if (SYSTEM.SDT[g].Drawing == 1):
            if (SYSTEM.SDT[g].Color == [0, 0, 0]):
                if (SYSTEM.SDT[g].Glass == 'MIRROR'):
                    LL_color = Mirror_color
                else:
                    LL_color = Glass_color
                if (SYSTEM.SDT[g].Glass == 'ABSORB'):
                    LL_color = Absorb_color
                color = LL_color
            else:
                color = SYSTEM.SDT[g].Color
            if (recorte == 1):
                clippedx = SYSTEM.BBB[n].clip('x', invert=False)
                p.add_mesh(clippedx, color ,opacity=OPA, smooth_shading=True, show_edges=None)
                clippedy = SYSTEM.BBB[n].clip('-y', invert=False)
                p.add_mesh(clippedy, color,opacity=OPA, smooth_shading=True, show_edges=None)
            if (recorte == 0):
                No_clipped = SYSTEM.BBB[n]
                p.add_mesh(No_clipped, color,opacity=OPA, smooth_shading=False, show_edges=None)
            if (recorte == 2):
                BBBB = SYSTEM.BBB[n]
                clippedx = BBBB.clip((1, 0, 0), invert=False)
                p.add_mesh(clippedx, color,opacity=OPA, smooth_shading=True, show_edges=None)
        n = (n + 1)

    if INST == True:
       for R in RAYS:
            if (len(R.RayWave) != 0):
                RW = np.asarray(R.RayWave)
                for i in range(0, np.shape(RW)[0]):
                    RGB = wavelength_to_rgb((RW[i] * 1000.0))
                    p.add_mesh(CCC[i], color=RGB, opacity=OPA, smooth_shading=True, line_width=1.0, show_edges=None)

    else:
        if (len(RAYS.RayWave) != 0):
            RW = np.asarray(RAYS.RayWave)
            for i in range(0, np.shape(RW)[0]):
                RGB = wavelength_to_rgb((RW[i] * 1000.0))
                p.add_mesh(CCC[i], color=RGB, opacity=OPA, smooth_shading=True, line_width=1.0, show_edges=None)

    p.add_axes(line_width=4)

    [cx,cy,cz]=p.center
    p.set_focus([cx,cy,cz])
    p.camera_position=[-1.0, 0.5, 1.0]

    p.set_viewup([0, 1.0, 0])
    p.enable_anti_aliasing()
    p.enable_eye_dome_lighting()
    p.disable_eye_dome_lighting()
    p.set_background(BackgCol, top=BackgColTop)

    p.add_text('KrakenOS',position="upper_left" ,font_size=20,color="royalblue")
    p.show_grid(font_size=6,color=GridCol)
    p.show(auto_close=False, interactive=True, interactive_update=False)

    [cpx,cpy,cpz]=p.camera_position

###############################################################################

def display3d(SYSTEM, RAYS, view=0, inline=False,     BackgCol= 'white', BackgColTop = 'white', GridCol="black", nrays = 0):

    """display3d.

    Parameters
    ----------
    SYSTEM :
        SYSTEM
    RAYS :
        RAYS
    view :
        view
    """

    ST1="KrakenOS v0.1. Executing Script: " + sys.argv[0]
    OPA = 0.95
    p = pv.Plotter(shape=(1, 1), title=ST1,notebook=inline)


    INST = isinstance(SYSTEM, list)
    if INST == False:
        SYSTEM = [SYSTEM]

    INST = isinstance(RAYS, list)
    if INST == False:
        RAYS = [RAYS]
    for system in SYSTEM:
        plot3d(system, view, p, OPA)
    for rays in RAYS:
        rayplot3d(rays, view, p, OPA, nrays)
    p.add_axes(line_width=4)

    [cx,cy,cz]=p.center
    p.set_focus([cx,cy,cz])
    p.camera_position=[-1.0, 0.5, 1.0]

    p.set_viewup([0, 1.0, 0])
    p.enable_anti_aliasing()
    p.enable_eye_dome_lighting()
    p.disable_eye_dome_lighting()
    p.set_background(BackgCol, top=BackgColTop)
    p.add_text('KrakenOS',position="upper_left" ,font_size=20,color="royalblue")
    p.show_grid(font_size=6,color=GridCol)
    p.show(auto_close=False, interactive=True, interactive_update=False)
    [cpx,cpy,cpz]=p.camera_position

###############################################################################

def rayplot3d(RAYS, view, p, OPA, nrays):
    CCC = pv.MultiBlock()
    for rays in RAYS.CC:
        RAY_VTK_OBJ = pv.lines_from_points(rays)
        CCC.append(RAY_VTK_OBJ)

    if (len(RAYS.RayWave) != 0):
        RW = np.asarray(RAYS.RayWave)
        if nrays == 0:
            NR = np.shape(RW)[0]
        else:
            NR = nrays

        for i in range(0, NR):
            RGB = wavelength_to_rgb((RW[i] * 1000.0))
            p.add_mesh(CCC[i], color=RGB, opacity=OPA, smooth_shading=True, line_width=1.0, show_edges=None)
    return 0

###############################################################################


def plot3d(SYSTEM, view, p, OPA):

    Absorb_color = np.array([(10 / 256), (23 / 256), (24 / 256)])
    Mirror_color = np.array([(189 / 256), (189 / 256), (189 / 256)])
    Glass_color=np.array([12/256, 238/256, 246/256])
    recorte = view
    NN = SYSTEM.AAA.n_blocks
    if (SYSTEM.SDT[0].Drawing == 0):
        points2 = np.c_[(0, 0, 0)]
        c = pv.PolyData(points2)
        cc = pv.PolyData(points2)
    if (SYSTEM.SDT[0].Drawing == 1):
        c = SYSTEM.AAA[0]
        cc = SYSTEM.AAA[0]

    for n in range(1, NN):
        if (SYSTEM.SDT[n].Drawing == 1):
            AAAA = SYSTEM.AAA[n]
            if (SYSTEM.SDT[n].Glass != 'NULL'):
                if (SYSTEM.SDT[n].Color == [0, 0, 0]):
                    if (SYSTEM.SDT[n].Glass == 'MIRROR'):
                        color = Mirror_color
                    else:
                        color = Glass_color
                    if (SYSTEM.SDT[n].Glass == 'ABSORB'):
                        color = Absorb_color
                else:
                    color = SYSTEM.SDT[n].Color

                cc = cc.merge(AAAA)
                if (recorte == 1):
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    clippedy = clippedx.clip((0, 1, 0), invert=False)
                    c = c.merge(clippedy)
                    clippedx = AAAA.clip(((- 1), 0, 0), invert=False)
                    clippedy = clippedx.clip((0, (- 1), 0), invert=False)
                    c = c.merge(clippedy)
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    clippedy = clippedx.clip((0, (- 1), 0), invert=False)
                    c = c.merge(clippedy)
                if (recorte == 0):
                    c = cc

                if (recorte == 2):
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    c = c.merge(clippedx)
                p.add_mesh(c, color, opacity=OPA, specular=1, specular_power=15, smooth_shading=True, show_edges=False)
                edges = c.extract_feature_edges(feature_angle=10, boundary_edges=True, feature_edges=False, manifold_edges=False)
                if SYSTEM.SDT[n].Solid_3d_stl != "None" and recorte ==0:
                    print(" ") # No edges
                else:
                    p.add_mesh(edges, 'red')
                points2 = np.c_[(0, 0, 0)]
                c = pv.PolyData(points2)
    p.add_mesh(SYSTEM.DDD, color=[0.5, 0.5, 0.5], opacity=OPA, show_edges=None)
    NN = SYSTEM.AAA.n_blocks
    n = 0
    for g in SYSTEM.side_number:
        if (SYSTEM.SDT[g].Drawing == 1):
            if (SYSTEM.SDT[g].Color == [0, 0, 0]):
                if (SYSTEM.SDT[g].Glass == 'MIRROR'):
                    LL_color = Mirror_color
                else:
                    LL_color = Glass_color
                if (SYSTEM.SDT[g].Glass == 'ABSORB'):
                    LL_color = Absorb_color
                color = LL_color
            else:
                color = SYSTEM.SDT[g].Color
            if (recorte == 1):
                clippedx = SYSTEM.BBB[n].clip('x', invert=False)
                p.add_mesh(clippedx, color ,opacity=OPA, smooth_shading=True, show_edges=None)
                clippedy = SYSTEM.BBB[n].clip('-y', invert=False)
                p.add_mesh(clippedy, color,opacity=OPA, smooth_shading=True, show_edges=None)
            if (recorte == 0):
                No_clipped = SYSTEM.BBB[n]
                p.add_mesh(No_clipped, color,opacity=OPA, smooth_shading=False, show_edges=None)
            if (recorte == 2):
                BBBB = SYSTEM.BBB[n]
                clippedx = BBBB.clip((1, 0, 0), invert=False)
                p.add_mesh(clippedx, color,opacity=OPA, smooth_shading=True, show_edges=None)
        n = (n + 1)
    return 0

###############################################################################
def display2d(SYSTEM, RAYS, view=0, arrow=0, nrays = 0):
    """display2d.

    Parameters
    ----------
    SYSTEM :
        SYSTEM
    RAYS :
        RAYS
    view :
        view
    """


    fig, ax1 = plt.subplots()

    INST = isinstance(SYSTEM, list)
    if INST == True:

        for SYS in SYSTEM:
            Plot2DSurf(SYS, view, ax1)
    else:
        Plot2DSurf(SYSTEM, view, ax1)

    INST = isinstance(RAYS, list)
    if INST == True:
        for R in RAYS:
            Plot2DRays(R, view, arrow, ax1, nrays)
    else:
        Plot2DRays(RAYS, view, arrow, ax1, nrays)

    plt.title('System Plot')
    plt.xlabel('Z (mm)')
    if (view == 0):
        plt.ylabel('Y (mm)')
    if (view == 1):
        plt.ylabel('X (mm)')
    plt.axis('equal')
    plt.show()

###############################################################################

def Plot2DSurf(SYSTEM, view, ax1):
    fs = 11
    NN = SYSTEM.AAA.n_blocks
    if (SYSTEM.SDT[0].Drawing == 0):
        points2 = np.c_[(0.0, 0.0, 0.0)]
        c = pv.PolyData(points2, force_float=False)
    if (SYSTEM.SDT[0].Drawing == 1):
        c = SYSTEM.AAA[0]
    sign=-1.0
    sn=0.65
    mx = []
    mn = []
    for n in range(0, NN):
        if SYSTEM.SDT[n].Solid_3d_stl != "None":
            solid = 1
        else:
            solid = 0

        sn = sn * sign
        if (SYSTEM.SDT[n].Drawing == 1):
            AAAA = SYSTEM.AAA[n]
            if (SYSTEM.SDT[n].Glass != 'NULL'):
                (PosX, PosY) = SYSTEM.SDT[n].Nm_Pos
                s = SYSTEM.SDT[n].Name
                ss = SYSTEM.Object_Num[n]
                if (view == 0):
                    LT = ''
                    (ax, ay, az) = edge_3d(AAAA, 1, 0, 0, solid)
                    (az, ay) = filter_face_2dplot(az, ay, solid)
                    ax1.plot(az, ay, LT, c='black', linewidth=0.5)
                    (ax, ay, az) = edge_3d(AAAA, (- 1), 0, 0, solid)
                    (az, ay) = filter_face_2dplot(az, ay, solid)
                    ax1.plot(az, ay, LT, c='black', linewidth=0.5)
                    ax1.text(((np.max(az) + PosX) + 1), ((np.max(ay) + PosY) - 1), s, fontsize=fs)
                    delta = ((np.max(ay) - np.min(ay)) / 10)


                    ABC=int(0.5*(np.argmin(ay) - np.argmax(ay)))
                    CordABC = az[ABC]


                    if SYSTEM.SDT[n].NumLabel == 1:

                        ax1.text(CordABC, sn*1.5*(np.min(ay) - (1.5 * delta)), (('[' + str(ss)) + ']'), fontsize=fs)
                        ax1.plot([CordABC, CordABC], [np.mean(ay), sn*1.5*(np.min(ay) - delta)], '-.', c='red', linewidth=0.5)
                    if ((PosX != 0) or (PosY != 0)):
                        ax1.arrow((np.max(az) + PosX), (np.max(ay) + PosY/2), (- PosX), (- PosY/2.0), head_width=0.5, head_length=1.0, fc='k', ec='k', length_includes_head=True)
                        ax1.arrow((np.max(az) + PosX), (np.max(ay) + PosY), (- PosX)*0, (- PosY)/2.0, head_width=0.1, head_length=0.0, fc='k', ec='k', length_includes_head=True)

                if (view == 1):
                    LT = ''
                    (ax, ay, az) = edge_3d(AAAA, 0, 1, 0, solid)
                    (az, ax) = filter_face_2dplot(az, ax, solid)
                    ax1.plot(az, ax, LT, c='black', linewidth=0.5)
                    (ax, ay, az) = edge_3d(AAAA, 0, (- 1), 0, solid)
                    (az, ax) = filter_face_2dplot(az, ax, solid)
                    ax1.plot(az, ax, LT, c='black', linewidth=0.5)
                    ax1.text(((np.max(az) + PosX) + 1), ((np.max(ax) + PosY) - 1), s, fontsize=fs)
                    delta = ((np.max(ax) - np.min(ax)) / 10)

                    ABC=int(0.5*(np.argmin(ax) - np.argmax(ax)))
                    CordABC = az[ABC]


                    if SYSTEM.SDT[n].NumLabel == 1:
                        ax1.text(CordABC, sn*1.5*(np.min(ax) - (1.5 * delta)), (('[' + str(ss)) + ']'), fontsize=fs)
                        ax1.plot([CordABC, CordABC], [np.mean(ax), sn*1.5*(np.min(ax) - delta)], '-.', c='red', linewidth=0.5)

                    if ((PosX != 0) or (PosY != 0)):
                        ax1.arrow((np.max(az) + PosX), (np.max(ax) + PosY/2), (- PosX), (- PosY/2.0), head_width=0.5, head_length=1.0, fc='k', ec='k', length_includes_head=True)
                        ax1.arrow((np.max(az) + PosX), (np.max(ax) + PosY), (- PosX)*0, (- PosY)/2.0, head_width=0.1, head_length=0.0, fc='k', ec='k', length_includes_head=True)

            mx.append(np.max(az))
            mn.append(np.min(az))

    NN = SYSTEM.BBB.n_blocks
    for n in range(0, NN):
        if SYSTEM.SDT[n].Solid_3d_stl != "None":
            solid = 1
        else:
            solid = 0

        TT = SYSTEM.BBB[n]


        sim = '-.'
        if (view == 0):
            (ax, ay, az) = edge_3d(TT, 1, 0, 0, solid)
            ax1.plot(az, ay, sim, c='black', linewidth=0.5)
            (ax, ay, az) = edge_3d(AAAA, (- 1), 0, 0, solid)
            ax1.plot(az, ay, sim, c='black', linewidth=0.5)
        if (view == 1):
            (ax, ay, az) = edge_3d(TT, 0, 1, 0, solid)
            ax1.plot(az, ax, sim, c='black', linewidth=0.5)
            (ax, ay, az) = edge_3d(AAAA, 0, (- 1), 0, solid)
            ax1.plot(az, ax, sim, c='black', linewidth=0.5)

    mx = np.asarray(mx)
    mn = np.asarray(mn)

    m1=np.min(mn)
    m2=np.max(mx)
    m3=(m2-m1)/20.0
    ax1.plot([m1-m3 , m2+m3],[0,0],'--', color="black", linewidth=1)

###############################################################################
def Plot2DRays(RAYS, view, arrow, ax1, nrays):

    CCC = pv.MultiBlock()
    for rays in RAYS.CC:
        RAY_VTK_OBJ = pv.lines_from_points(rays)
        CCC.append(RAY_VTK_OBJ)

    if (len(RAYS.RayWave) != 0):
        RW = np.asarray(RAYS.RayWave)
        if nrays == 0:
            NR = np.shape(RW)[0]
        else:
            NR = nrays
        for i in range(0, NR):
            RGB = wavelength_to_rgb((RW[i] * 1000.0))
            RRR = CCC[i]
            Ax = RRR.points[:, 0]
            Ay = RRR.points[:, 1]
            Az = RRR.points[:, 2]
            for f in range(1,len(Ax)):
                if (view == 0):
                    line = ax1.plot([Az[f-1],Az[f]], [Ay[f-1],Ay[f]], color=RGB, linewidth=0.5)[0]
                if (view == 1):
                    line = ax1.plot([Az[f-1],Az[f]], [Ax[f-1],Ax[f]], color=RGB, linewidth=0.5)[0]
                if arrow != 0:
                    add_arrow(line, size=5*arrow)

###############################################################################

def edge_3d(MeshObject, cx, cy, xz, solid):
    """edge_3d.

    Parameters
    ----------
    MeshObject :
        MeshObject
    cx :
        cx
    cy :
        cy
    xz :
        xz
    """
    c = MeshObject.clip((cx, cy, xz), invert=False)
    edges = c.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False)

    if solid ==0:
        Ax = edges.points[:, 0]
        Ay = edges.points[:, 1]
        Az = edges.points[:, 2]

    else:
        Ax = c.points[:, 0]
        Ay = c.points[:, 1]
        Az = c.points[:, 2]

    Ax = np.asarray(Ax)
    Ay = np.asarray(Ay)
    Az = np.asarray(Az)
    Xe = []
    Ye = []
    Ze = []
    if (cx == 1):
        i = np.argmin(Ax)
    if (cx == (- 1)):
        i = np.argmax(Ax)
    if (cy == 1):
        i = np.argmin(Ay)
    if (cy == (- 1)):
        i = np.argmax(Ay)
    x0 = Ax[i]
    y0 = Ay[i]
    z0 = Az[i]
    for j in range(0, (np.shape(Ax)[0] - 1)):
        AAx = Ax[i]
        AAy = Ay[i]
        AAz = Az[i]
        Xe.append(AAx)
        Ye.append(AAy)
        Ze.append(AAz)
        Ax = np.delete(Ax, i)
        Ay = np.delete(Ay, i)
        Az = np.delete(Az, i)
        X = (Ax - AAx)
        Y = (Ay - AAy)
        Z = (Az - AAz)
        R = np.sqrt((((X * X) + (Y * Y)) + (Z * Z)))
        i = np.argmin(R)
    Xe.append(x0)
    Ye.append(y0)
    Ze.append(z0)
    Xe = np.asarray(Xe)
    Ye = np.asarray(Ye)
    Ze = np.asarray(Ze)
    return (Xe, Ye, Ze)

def filter_face_2dplot(v1, v2, solid):
    """filter_face_2dplot.

    Parameters
    ----------
    v1 :
        v1
    v2 :
        v2
    """
    if solid == 0:
        av1 = np.copy(v1)
        av2 = np.copy(v2)
        av1 = np.roll(av1, (- 1))
        av2 = np.roll(av2, (- 1))
        yy = (v1 - av1)
        zz = (v2 - av2)
        R = np.sqrt(((yy * yy) + (zz * zz)))
        M = np.mean(R)
        AW = np.argwhere((R < M))
        v1 = v1[AW]
        v2 = v2[AW]
        av1 = np.copy(v1)
        av2 = np.copy(v2)
        av1 = np.roll(av1, 1)
        av2 = np.roll(av2, 1)
        yy = (v1 - av1)
        zz = (v2 - av2)
        R = np.sqrt(((yy * yy) + (zz * zz)))
        AW = np.argmax(R)
        v1 = np.roll(v1, (- AW))
        v2 = np.roll(v2, (- AW))

    return (v1, v2)

#####################################

class display3d_OB():
    def __init__(self):


        self.view=0
        self.inline=False
        self.BackgCol= 'white'
        self.BackgColTop = 'white'
        self.GridCol="black"
        self.ST1 = " "
        self.p = pv.Plotter(shape=(1, 1), title=self.ST1,notebook = self.inline)
        self.SYSTEM = []
        self.RAYS = []

    def plot(self):
        display3d_4OB(self.SYSTEM, self.RAYS, self.view, self.inline, self.BackgCol, self.BackgColTop, self.GridCol, self.p)

def display3d_4OB(SYSTEM, RAYS, view, inline, BackgCol, BackgColTop, GridCol, p):


    """display3d.

    Parameters
    ----------
    SYSTEM :
        SYSTEM
    RAYS :
        RAYS
    view :
        view
    """

    OPA = 0.95

    CCC = pv.MultiBlock()
    INST = isinstance(RAYS, list)
    if INST == True:
        for R in RAYS:
            for rays in R.CC:
                RAY_VTK_OBJ = pv.lines_from_points(rays)
                CCC.append(RAY_VTK_OBJ)
    else:
        for rays in RAYS.CC:
            RAY_VTK_OBJ = pv.lines_from_points(rays)
            CCC.append(RAY_VTK_OBJ)

    Absorb_color = np.array([(10 / 256), (23 / 256), (24 / 256)])
    Mirror_color = np.array([(189 / 256), (189 / 256), (189 / 256)])
    Glass_color=np.array([12/256, 238/256, 246/256])
    recorte = view
    NN = SYSTEM.AAA.n_blocks
    if (SYSTEM.SDT[0].Drawing == 0):
        points2 = np.c_[(0, 0, 0)]
        c = pv.PolyData(points2)
        cc = pv.PolyData(points2)
    if (SYSTEM.SDT[0].Drawing == 1):
        c = SYSTEM.AAA[0]
        cc = SYSTEM.AAA[0]

    for n in range(1, NN):
        if (SYSTEM.SDT[n].Drawing == 1):
            AAAA = SYSTEM.AAA[n]
            if (SYSTEM.SDT[n].Glass != 'NULL'):
                if (SYSTEM.SDT[n].Color == [0, 0, 0]):
                    if (SYSTEM.SDT[n].Glass == 'MIRROR'):
                        color = Mirror_color
                    else:
                        color = Glass_color
                    if (SYSTEM.SDT[n].Glass == 'ABSORB'):
                        color = Absorb_color
                else:
                    color = SYSTEM.SDT[n].Color

                cc = cc.merge(AAAA)
                if (recorte == 1):
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    clippedy = clippedx.clip((0, 1, 0), invert=False)
                    c = c.merge(clippedy)
                    clippedx = AAAA.clip(((- 1), 0, 0), invert=False)
                    clippedy = clippedx.clip((0, (- 1), 0), invert=False)
                    c = c.merge(clippedy)
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    clippedy = clippedx.clip((0, (- 1), 0), invert=False)
                    c = c.merge(clippedy)
                if (recorte == 0):
                    c = cc

                if (recorte == 2):
                    clippedx = AAAA.clip((1, 0, 0), invert=False)
                    c = c.merge(clippedx)
                p.add_mesh(c, color, opacity=OPA, specular=1, specular_power=15, smooth_shading=True, show_edges=False)
                edges = c.extract_feature_edges(feature_angle=10, boundary_edges=True, feature_edges=False, manifold_edges=False)
                if SYSTEM.SDT[n].Solid_3d_stl != "None" and recorte ==0:
                    print(" ") # No edges
                else:
                    p.add_mesh(edges, 'red')
                points2 = np.c_[(0, 0, 0)]
                c = pv.PolyData(points2)
    p.add_mesh(SYSTEM.DDD, color=[0.5, 0.5, 0.5], opacity=OPA, show_edges=None)
    NN = SYSTEM.AAA.n_blocks
    n = 0
    for g in SYSTEM.side_number:
        if (SYSTEM.SDT[g].Drawing == 1):
            if (SYSTEM.SDT[g].Color == [0, 0, 0]):
                if (SYSTEM.SDT[g].Glass == 'MIRROR'):
                    LL_color = Mirror_color
                else:
                    LL_color = Glass_color
                if (SYSTEM.SDT[g].Glass == 'ABSORB'):
                    LL_color = Absorb_color
                color = LL_color
            else:
                color = SYSTEM.SDT[g].Color
            if (recorte == 1):
                clippedx = SYSTEM.BBB[n].clip('x', invert=False)
                p.add_mesh(clippedx, color ,opacity=OPA, smooth_shading=True, show_edges=None)
                clippedy = SYSTEM.BBB[n].clip('-y', invert=False)
                p.add_mesh(clippedy, color,opacity=OPA, smooth_shading=True, show_edges=None)
            if (recorte == 0):
                No_clipped = SYSTEM.BBB[n]
                p.add_mesh(No_clipped, color,opacity=OPA, smooth_shading=False, show_edges=None)
            if (recorte == 2):
                BBBB = SYSTEM.BBB[n]
                clippedx = BBBB.clip((1, 0, 0), invert=False)
                p.add_mesh(clippedx, color,opacity=OPA, smooth_shading=True, show_edges=None)
        n = (n + 1)

    if INST == True:
       for R in RAYS:
            if (len(R.RayWave) != 0):
                RW = np.asarray(R.RayWave)
                for i in range(0, np.shape(RW)[0]):
                    RGB = wavelength_to_rgb((RW[i] * 1000.0))
                    p.add_mesh(CCC[i], color=RGB, opacity=OPA, smooth_shading=True, line_width=1.0, show_edges=None)

    else:
        if (len(RAYS.RayWave) != 0):
            RW = np.asarray(RAYS.RayWave)
            for i in range(0, np.shape(RW)[0]):
                RGB = wavelength_to_rgb((RW[i] * 1000.0))
                p.add_mesh(CCC[i], color=RGB, opacity=OPA, smooth_shading=True, line_width=1.0, show_edges=None)

    p.add_axes(line_width=4)
    [cx,cy,cz]=p.center
    p.set_focus([cx,cy,cz])
    p.camera_position=[-1.0, 0.5, 1.0]
    p.set_viewup([0, 1.0, 0])
    p.enable_anti_aliasing()
    p.enable_eye_dome_lighting()
    p.disable_eye_dome_lighting()
    p.set_background(BackgCol, top=BackgColTop)
    p.add_text('KrakenOS',position="upper_left" ,font_size=20,color="royalblue")
    p.show_grid(font_size=6,color=GridCol)
    p.show(auto_close=False, interactive=True, interactive_update=True)
    [cpx,cpy,cpz]=p.camera_position
