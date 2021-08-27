# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 17:42:47 2021

@author: JOELHERRERAVAZQUEZ
"""

import csv
import numpy as np


##############################################################################

def load_alluminum_complex(file):
    """ hsahdasjhgd """

    w = []
    n = []
    k = []
    with open(file, 'rt') as f:
        csv_reader = csv.reader(f, delimiter=';')
        next(csv_reader)  # skip the heading
        for line in csv_reader:
            w.append(float(line[0]))
            n.append(float(line[1]))
            k.append(float(line[2]))

    return w, n, k


# ##############################################################################

###########################################################################
# def load_Catalog2(FileCat):
#     global NM, ED, LD, OD, TD, CD, file
#     ARR_CAT = []
#     ARR_NM = []
#     ARR_ED = []
#     ARR_CD = []
#     ARR_TD = []
#     ARR_OD = []
#     ARR_LD = []
#     ARR_IT = []
#     cont = 0
#     names = []
#     ran = []
#     cat = []

#     for file in FileCat:
#         print(file)
#         f = open(file, "r")
#         NM_FileValidity = []
#         for x in f:
#             cat.append(x)
#             st = x[0:2]

#             if st == "NM":
#                 # print(cont,st)
#                 spt = x.split()
#                 names.append(spt[1])
#                 NM_FileValidity.append(spt[1])
#                 ran.append(float(cont))
#             cont = cont + 1
#         ran.append(float(cont))

#         valid = len(NM_FileValidity)
#         if valid == 0:
#             Message = "The glass catalog (" + file + ") is not valid or is not utf8 encoding"
#             raise ValueError(Message)

#     names = np.asarray(names)
#     ran = np.asarray(ran)

#     for i in range(0, np.shape(names)[0]):
#         Glass = names[i]

#         place = np.argwhere(names == Glass)
#         ITT = []
#         for j in range(int(ran[place][0][0]), int(ran[place + 1][0][0])):
#             cadena = cat[j].split()
#             cad = cat[j][2:].split()

#             if cadena[0] == "NM":
#                 NM = cad
#             if cadena[0] == "ED":
#                 ED = cad
#             if cadena[0] == "CD":
#                 CD = cad
#             if cadena[0] == "TD":
#                 TD = cad
#             if cadena[0] == "OD":
#                 OD = cad
#             if cadena[0] == "LD":
#                 LD = cad
#             if cadena[0] == "IT":
#                 IT = cad
#                 ITT.append(IT)

#         NM = np.asarray(NM[1:-1], dtype=np.float64)
#         ED = np.asarray(ED, dtype=np.float64)
#         CD = np.asarray(CD, dtype=np.float64)
#         TD = np.asarray(TD, dtype=np.float64)
#         OD = np.asarray(OD)
#         LD = np.asarray(LD, dtype=np.float64)
#         IT = np.asarray(ITT, dtype=np.float64).T

#         ARR_CAT.append(file)
#         ARR_NM.append(NM)
#         ARR_ED.append(ED)
#         ARR_CD.append(CD)
#         ARR_TD.append(TD)
#         ARR_OD.append(OD)
#         ARR_LD.append(LD)
#         ARR_IT.append(IT)

#     CATALOG = [ARR_CAT, names, ARR_NM, ARR_ED, ARR_CD, ARR_TD, ARR_OD, ARR_LD, ARR_IT]
#     return CATALOG

# ##############################################################################


def load_Catalog(FileCat):
    cat = []
    ARR_CAT = []
    ARR_NM = []
    ARR_ED = []
    ARR_CD = []
    ARR_TD = []
    ARR_OD = []
    ARR_LD = []
    ARR_IT = []
    print("Loading glass calatogs:")
    for file in FileCat:
        ARR_CAT.append(file)
        f = open(file, "r")
        for x in f:
            cat.append(x)

    con = 0
    coords = []
    names = []
    for f in cat:
        cadena = f.split()
        cad = cadena[0]
        if cad == "NM":
            coords.append(con)
            names.append(cadena[1])
            # print(cad, con, cadena[1])
        con = con + 1

    names = np.asarray(names)

    for i in range(0, len(coords) - 1):
        # print(names[i], coords[i], coords[i+1])

        ITT = []
        for j in range(coords[i], coords[i + 1]):
            cadena = cat[j].split()
            cad = cat[j][2:].split()

            if cadena[0] == "NM":
                NM = cad
            if cadena[0] == "ED":
                ED = cad
            if cadena[0] == "CD":
                CD = cad
            if cadena[0] == "TD":
                TD = cad
            if cadena[0] == "OD":
                OD = cad
            if cadena[0] == "LD":
                LD = cad
            if cadena[0] == "IT":
                IT = cad
                ITT.append(IT)

        NM = np.asarray(NM[1:-1], dtype=np.float64)
        ED = np.asarray(ED, dtype=np.float64)
        CD = np.asarray(CD, dtype=np.float64)
        TD = np.asarray(TD, dtype=np.float64)
        OD = np.asarray(OD)
        LD = np.asarray(LD, dtype=np.float64)
        IT = np.asarray(ITT, dtype=np.float64).T

        ARR_NM.append(NM)
        ARR_ED.append(ED)
        ARR_CD.append(CD)
        ARR_TD.append(TD)
        ARR_OD.append(OD)
        ARR_LD.append(LD)
        ARR_IT.append(IT)

    CATALOG = [ARR_CAT, names, ARR_NM, ARR_ED, ARR_CD, ARR_TD, ARR_OD, ARR_LD, ARR_IT]
    return CATALOG
