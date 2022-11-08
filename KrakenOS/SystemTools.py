import csv
import numpy as np

def load_alluminum_complex(file):
    """load_alluminum_complex.

    Parameters
    ----------
    file :
        file
    """
    w = []
    n = []
    k = []
    with open(file, 'rt') as f:
        csv_reader = csv.reader(f, delimiter=';')
        next(csv_reader)
        for line in csv_reader:
            w.append(float(line[0]))
            n.append(float(line[1]))
            k.append(float(line[2]))
    return (w, n, k)

def load_Catalog(FileCat):
    """load_Catalog.

    Parameters
    ----------
    FileCat :
        FileCat
    """
    cat = []
    ARR_CAT = []
    ARR_NM = []
    ARR_ED = []
    ARR_CD = []
    ARR_TD = []
    ARR_OD = []
    ARR_LD = []
    ARR_IT = []
    print('Processing glass catalog:')
    for file in FileCat:
        ARR_CAT.append(file)
        f = open(file, 'r')
        for x in f:

            if not x.isspace():
                 cat.append(x)
    con = 0
    coords = []
    names = []
    for f in cat:
        cadena = f.split()
        cad = cadena[0]
        if (cad == 'NM'):
            coords.append(con)
            names.append(cadena[1])
        con = (con + 1)
    names = np.asarray(names)
    for i in range(0, (len(coords) - 1)):
        ITT = []
        for j in range(coords[i], coords[(i + 1)]):
            cadena = cat[j].split()
            cad = cat[j][2:].split()
            if (cadena[0] == 'NM'):
                NM = cad
            if (cadena[0] == 'ED'):
                ED = cad
                if ED[1]=="-":
                    ED[1]="0.0"
            if (cadena[0] == 'CD'):
                CD = cad
            if (cadena[0] == 'TD'):
                TD = cad
            if (cadena[0] == 'OD'):
                OD = cad
            if (cadena[0] == 'LD'):
                LD = cad
            if (cadena[0] == 'IT'):
                IT = cad
                if len(IT)==3:
                    ITT.append(IT)

        NM = np.asarray(NM[1:(- 1)], dtype=np.float64)
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
