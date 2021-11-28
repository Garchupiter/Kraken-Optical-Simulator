
import numpy as np
import KrakenOS as Kos
import scipy
import copy













def Phase(PUPIL):
    """Phase.

    Parameters
    ----------
    PUPIL :
        PUPIL
    """


#############################################################
    """ Crea una pupila con las propiedades pasadas como parámetro"""


    # alpha = abcd[0]
    # betha = abcd[1]
    # gamma = abcd[2]
    # ro = abcd[3]


    SYSTEM = PUPIL.SYSTEM
    sup = PUPIL.Surf
    W = PUPIL.W
    ApType = PUPIL.ApertureType
    ApVal = PUPIL.ApertureValue
    configuracion_1 = SYSTEM.SETUP
    W = W

    """ Se calcula la pupila """
    Pup = Kos.PupilCalc(SYSTEM, sup, W, ApType, ApVal)

    [PIX, PIY, PIZ] = PUPIL.PosPupInp
    [POX, POY, POZ] = PUPIL.PosPupOutFoc

    [PIL, PIM, PIN] = PUPIL.DirPupSal
    [POL, POM, PON] = PUPIL.DirPupSal

    Pup.Samp = PUPIL.Samp
    Pup.Ptype = PUPIL.Ptype
    Pup.FieldY = PUPIL.FieldY
    Pup.FieldX = PUPIL.FieldX
    Pup.FieldType = PUPIL.FieldType

    """se generan los rayos para trazar por la pupila"""
    x, y, z, L, M, N = Pup.Pattern2Field()





    # Pup.Samp = 1
    Pup.Ptype = 'chief'
    xc, yc, zc, Lc, Mc, Nc = Pup.Pattern2Field()

    pSource_0_c = [xc[0], yc[0], zc[0]]
    dCos_c = [Lc[0], Mc[0], Nc[0]]

    """ Traza el rayo principal"""
    SYSTEM.IgnoreVignetting()
    SYSTEM.Trace(pSource_0_c, dCos_c, W)

    """ Toma las coordenadas y los cosenos directores del rayo
    principal a la salida del sistema"""

    [Xcc, Ycc, Zcc] = SYSTEM.XYZ[(- 1)]
    [Lcc, Mcc, Ncc] = SYSTEM.LMN[(- 1)]


    resf_out =  np.sqrt(((POX - Xcc)**2.0) + ((POY- Ycc)**2.0) + (POZ**2.0))

# - - - - - - - - - - - - - - - - - - - - - - - - -
    ArrSup = PUPIL.SYSTEM.SDT

    SS = []
    for ii in range(0,len(ArrSup)):
        ArrSup[ii].EraseVTK()
        ss = copy.deepcopy(ArrSup[ii])
        SS.append(ss)

    DisPup = SS[0].Thickness

# - - - - - - - - - - - - - - - - - - - - - -
    """ No mover """


    """ No mover """
    INP_P = Kos.surf()
    INP_P.Rc = 0
    INP_P.Thickness = DisPup
    INP_P.Diameter = PUPIL.RadPupInp * 3.0
    INP_P.Glass = 'AIR'
    INP_P.DespX = xc[0]
    INP_P.DespY = yc[0]
    INP_P.TiltX = np.rad2deg(np.arcsin((- Mc[0])))
    INP_P.TiltY = np.rad2deg(np.arcsin((Lc[0] / np.cos(np.arcsin((- Mc[0]))))))
    INP_P.AxisMove = 0
    INP_P.Order = 1

# - - - - - - - - - - - - - - - - - - - - - -

    OUT_P = Kos.surf()
    OUT_P.Rc = -np.sign(POZ)*np.abs(resf_out) #+ro
    OUT_P.Thickness = 0
    OUT_P.Diameter = PUPIL.RadPupOut * 3
    OUT_P.Glass = 'AIR'
    OUT_P.DespX = POX #+ alpha
    OUT_P.DespY = POY #+ betha
    OUT_P.AxisMove = 0
    OUT_P.Order = 1

    OUT_P.TiltX = np.rad2deg(np.arcsin((- Mcc)))
    OUT_P.TiltY = np.rad2deg(np.arcsin((Lcc / np.cos(np.arcsin((- Mcc))))))

    Test = Kos.surf()
    Test.Diameter = PUPIL.RadPupOut * 3

# - - - - - - - - - - - - - - - - - - - - - - - - -

    SS.insert(1, INP_P)
    SS.append(OUT_P)


    SS[0].Thickness = 0.0
    SS[-2].Thickness = POZ #+ gamma


    SYSTEM = Kos.system(SS, configuracion_1)

# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Se configuran los contenedores de rayos"""
    RR = Kos.raykeeper(SYSTEM)

# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Configura al sistema para ignorar los viñeteos """
    SYSTEM.IgnoreVignetting()

# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Traza los rayos por systema"""

    for i in range(0, len(x)):
        """ Traza los rayos al diametro de la pupila"""
        pSource_0 = [x[i], y[i], z[i]]
        dCos = [L[i], M[i], N[i]]
        SYSTEM.Trace(pSource_0, dCos, W)
        RR.push()

# - - - - - - - - - - - - - - - - - - - - - - - - -

    All_Rays = np.asarray(RR.OP)

    axis = 1


    AR = np.sum(All_Rays,axis) - All_Rays[:,0] - (All_Rays[:,-1]*2.)

# - - - - - - - - - - - - - - - - - - - - - - - - -

    """ Traza el rayo principal"""

    SYSTEM.Trace(pSource_0_c, dCos_c, W)

    Chief = np.asarray(SYSTEM.OP)

    CH = np.sum(Chief) - Chief[0] - (Chief[-1]*2.)

# - - - - - - - - - - - - - - - - - - - - - - - - -

    Wi = (CH - AR )


# - - - - - - - - - - - - - - - - - - - - - - - - -

    XPUP = x - xc[0]
    YPUP = y - yc[0]


    x0, y0, z0, k0, l0, mi = RR.pick(1)
    xi, yi, zi, ki, li, mi = RR.pick(-1)



    # V = np.sum((Wi - (ki * alpha) - (li * betha) - gamma - ((mi - 1) * ro))**2.0)/ (len(Wi))

    # A1 = np.sum(ki ** 2.0)
    # A2 = np.sum(ki * li)
    # A3 = np.sum(ki)
    # A4 = np.sum(ki * (mi - 1.0))

    # B1 = np.sum(li * ki)
    # B2 = np.sum(li * li)
    # B3 = np.sum(li)
    # B4 = np.sum(li * (mi - 1.0))

    # C1 = np.sum(ki)
    # C2 = np.sum(li)
    # C3 = len(xi)
    # C4 = np.sum((mi - 1.0))

    # D1 = np.sum((mi - 1.0) * ki)
    # D2 = np.sum((mi - 1.0) * li)
    # D3 = np.sum((mi - 1.0))
    # D4 = np.sum((mi - 1.0) ** 2.0)



    # d1 = np.sum(Wi * ki)
    # d2 = np.sum(Wi * li)
    # d3 = np.sum(Wi)
    # d4 = np.sum(Wi * (mi - 1.0))


    # A = np.array([
    #     [A1, B1, C1, D1],
    #     [A2, B2, C2, D2],
    #     [A3, B3, C3, D3],
    #     [A4, B4, C4, D4]
    #     ])


    # d = np.array([d1, d2, d3, d4])

    # abcd = np.linalg.solve(A, d)


    # alpha = abcd[0]
    # betha = abcd[1]
    # gamma = abcd[2]
    # ro = abcd[3]


# -------------------------------------------------

    Wi = ((Wi * 1000.0) / W)

    P2V = (np.max(Wi) - np.min(Wi))


    RR.clean()
    RR.push()
    # Kos.display2d(SYSTEM,RR,1,0)



    return ((YPUP / Pup.RadPupInp), (XPUP / Pup.RadPupInp), Wi, np.abs(P2V))








# def Phase2(PUPIL):

#     """Phase2.
#     in development

#     Parameters
#     ----------
#     PUPIL :
#         PUPIL
#     """


# #############################################################
#     """ Crea una pupila con las propiedades pasadas como parámetro"""


#     # alpha = abcd[0]
#     # betha = abcd[1]
#     # gamma = abcd[2]
#     # ro = abcd[3]


#     SYSTEM = PUPIL.SYSTEM
#     sup = PUPIL.Surf
#     W = PUPIL.W
#     ApType = PUPIL.ApertureType
#     ApVal = PUPIL.ApertureValue
#     configuracion_1 = SYSTEM.SETUP
#     W = W

#     """ Se calcula la pupila """
#     Pup = Kos.PupilCalc(SYSTEM, sup, W, ApType, ApVal)

#     [PIX, PIY, PIZ] = PUPIL.PosPupInp
#     [POX, POY, POZ] = PUPIL.PosPupOutFoc

#     [PIL, PIM, PIN] = PUPIL.DirPupSal
#     [POL, POM, PON] = PUPIL.DirPupSal

#     Pup.Samp = PUPIL.Samp
#     Pup.Ptype = PUPIL.Ptype
#     Pup.FieldY = PUPIL.FieldY
#     Pup.FieldX = PUPIL.FieldX
#     Pup.FieldType = PUPIL.FieldType

#     """se generan los rayos para trazar por la pupila"""
#     x, y, z, L, M, N = Pup.Pattern2Field()


#     # Pup.Samp = 1
#     Pup.Ptype = 'chief'
#     xc, yc, zc, Lc, Mc, Nc = Pup.Pattern2Field()

#     pSource_0_c = [xc[0], yc[0], zc[0]]
#     dCos_c = [Lc[0], Mc[0], Nc[0]]

#     """ Traza el rayo principal"""
#     SYSTEM.IgnoreVignetting()
#     SYSTEM.Trace(pSource_0_c, dCos_c, W)

#     """ Toma las coordenadas y los cosenos directores del rayo
#     principal a la salida del sistema"""

#     [Xcc, Ycc, Zcc] = SYSTEM.XYZ[(- 1)]
#     [Lcc, Mcc, Ncc] = SYSTEM.LMN[(- 1)]


#     resf_out =  np.sqrt(((POX - Xcc)**2.0) + ((POY- Ycc)**2.0) + (POZ**2.0))

# # - - - - - - - - - - - - - - - - - - - - - - - - -


#     ArrSup = PUPIL.SYSTEM.SDT

#     SS = []
#     for ii in range(0,len(ArrSup)):
#         ArrSup[ii].EraseVTK()
#         ss = copy.deepcopy(ArrSup[ii])
#         SS.append(ss)

#     DisPup = SS[0].Thickness

# # - - - - - - - - - - - - - - - - - - - - - -
#     """ No mover """
#     SS[0].Thickness = 0.0
#     SS[-1].Thickness = np.sign(POZ)*np.abs(resf_out) #+ gamma

# # - - - - - - - - - - - - - - - - - - - - - -

#     """ No mover """
#     INP_P = Kos.surf()
#     INP_P.Rc = 0
#     INP_P.Thickness = DisPup
#     INP_P.Diameter = PUPIL.RadPupInp * 3.0
#     INP_P.Glass = 'AIR'
#     INP_P.DespX = xc[0]
#     INP_P.DespY = yc[0]
#     INP_P.TiltX = np.rad2deg(np.arcsin((- Mc[0])))
#     INP_P.TiltY = np.rad2deg(np.arcsin((Lc[0] / np.cos(np.arcsin((- Mc[0]))))))
#     INP_P.AxisMove = 0
#     INP_P.Order = 1

# # - - - - - - - - - - - - - - - - - - - - - -

#     OUT_P = Kos.surf()
#     OUT_P.Rc = -np.sign(POZ)*np.abs(resf_out) #+ ro
#     OUT_P.Thickness = 0
#     OUT_P.Diameter = PUPIL.RadPupOut * 3
#     OUT_P.Glass = 'AIR'
#     OUT_P.DespX = Xcc #+ alpha
#     OUT_P.DespY = Ycc #+ betha
#     OUT_P.AxisMove = 0
#     OUT_P.Order = 1


#     Test = Kos.surf()
#     Test.Diameter = PUPIL.RadPupOut * 3

# # - - - - - - - - - - - - - - - - - - - - - - - - -


#     SS.insert(1, INP_P)
#     SS.append(OUT_P)
#     SYSTEM = Kos.system(SS, configuracion_1)

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Se configuran los contenedores de rayos"""
#     RR = Kos.raykeeper(SYSTEM)

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Configura al sistema para ignorar los viñeteos """
#     SYSTEM.IgnoreVignetting()

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Traza los rayos por systema"""

#     for i in range(0, len(x)):
#         """ Traza los rayos al diametro de la pupila"""
#         pSource_0 = [x[i], y[i], z[i]]
#         dCos = [L[i], M[i], N[i]]
#         SYSTEM.Trace(pSource_0, dCos, W)
#         RR.push()

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     All_Rays = np.asarray(RR.OP)

#     axis = 1
#     AR = (np.sum(All_Rays,axis) - All_Rays[:,0]) - All_Rays[:,-1]*2


# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     """ Traza el rayo principal"""

#     SYSTEM.Trace(pSource_0_c, dCos_c, W)

#     Chief = np.asarray(SYSTEM.OP)

#     CH = (np.sum(Chief) - Chief[0]) - Chief[-1]*2

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     Wi = CH - AR

# # - - - - - - - - - - - - - - - - - - - - - - - - -

#     XPUP = x - xc[0]
#     YPUP = y - yc[0]




# # -------------------------------------------------



#     x0, y0, z0, k0, l0, mi = RR.pick(1)
#     xi, yi, zi, ki, li, mi = RR.pick(-1)
#     # V = np.sum((Wi - (ki * alpha) - (li * betha) - gamma - ((mi - 1) * ro))**2.0)/ (len(Wi))

#     # A1 = np.sum(ki ** 2.0)
#     # A2 = np.sum(ki * li)
#     # A3 = np.sum(ki)
#     # A4 = np.sum(ki * (mi - 1.0))

#     # B1 = np.sum(li * ki)
#     # B2 = np.sum(li * li)
#     # B3 = np.sum(li)
#     # B4 = np.sum(li * (mi - 1.0))

#     # C1 = np.sum(ki)
#     # C2 = np.sum(li)
#     # C3 = len(xi)
#     # C4 = np.sum((mi - 1.0))

#     # D1 = np.sum((mi - 1.0) * ki)
#     # D2 = np.sum((mi - 1.0) * li)
#     # D3 = np.sum((mi - 1.0))
#     # D4 = np.sum((mi - 1.0) ** 2.0)



#     # d1 = np.sum(Wi * ki)
#     # d2 = np.sum(Wi * li)
#     # d3 = np.sum(Wi)
#     # d4 = np.sum(Wi * (mi - 1.0))


#     # A = np.array([
#     #     [A1, B1, C1, D1],
#     #     [A2, B2, C2, D2],
#     #     [A3, B3, C3, D3],
#     #     [A4, B4, C4, D4]
#     #     ])


#     # d = np.array([d1, d2, d3, d4])

#     # abcd = np.linalg.solve(A, d)





#     # alpha = abcd[0]
#     # betha = abcd[1]
#     # gamma = abcd[2]
#     # ro = abcd[3]



# # -------------------------------------------------

#     Wi = ((Wi * 1000.0) / W)




#     P2V = (np.max(Wi) - np.min(Wi))

#     RR.clean()
#     RR.push()
#     # Kos.display2d(SYSTEM,RR,1,0)


#     # Wi = Wi - np.median(Wi)

#     return ((YPUP / Pup.RadPupInp), (XPUP / Pup.RadPupInp), Wi, np.abs(P2V))


