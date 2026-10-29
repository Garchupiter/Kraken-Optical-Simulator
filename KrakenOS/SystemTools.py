import csv
import numpy as np



def _looks_numeric(token):
    try:
        float(token)
    except Exception:
        return False
    return True


def load_metal_complex(file):
    """load_metal_complex.

    Parameters
    ----------
    file :
        file
    """
    w = []
    n = []
    k = []
    with open(file, 'rt', encoding='utf-8') as f:
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

    print('Loading glass calatogs:')
    for file in FileCat:
        ARR_CAT.append(file)
        try:
            f = open(file, 'r', encoding='UTF8')
            for x in f:
                if not x.isspace():
                     cat.append(x)
        except:
            f = open(file, 'r', encoding='UTF16')
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
    section_keys = {"NM", "ED", "CD", "TD", "OD", "LD", "IT"}
    continuation_keys = {"ED", "CD", "TD", "OD", "LD"}
    # each glass's block runs to the next NM line, the last one to the end of the catalogs
    coords.append(len(cat))
    for i in range(0, (len(coords) - 1)):
        
        ITT = []
        NM = []
        ED = []
        CD = []
        TD = []
        OD = []
        LD = []
        IT = []
        current_key = None
        for j in range(coords[i], coords[(i + 1)]):
            cadena = cat[j].split()
            if not cadena:
                continue
            # a section's data may wrap onto following lines (Nikon's CD and TD do): a line that
            # starts with a number continues the section above it
            if cadena[0] in section_keys:
                current_key = cadena[0]
                cad = cadena[1:]
            else:
                if current_key not in continuation_keys or not _looks_numeric(cadena[0]):
                    continue
                cad = cadena
            if current_key == 'NM' and cadena[0] == 'NM':
                NM = cad
            if current_key == 'ED':
                ED = cad if cadena[0] == 'ED' else ED + cad
                if len(ED) > 1 and ED[1] == "-":
                    ED[1] = "0.0"
            if current_key == 'CD':
                CD = cad if cadena[0] == 'CD' else CD + cad
            if current_key == 'TD':
                TD = cad if cadena[0] == 'TD' else TD + cad
            if current_key == 'OD':
                OD = cad if cadena[0] == 'OD' else OD + cad
            if current_key == 'LD':
                LD = cad if cadena[0] == 'LD' else LD + cad
            if current_key == 'IT' and cadena[0] == 'IT':
                IT = cad
                if len(IT) == 3:
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
