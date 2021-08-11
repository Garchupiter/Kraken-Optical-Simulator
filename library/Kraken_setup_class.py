# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 17:55:21 2021

@author: JOELHERRERAVAZQUEZ
"""

import os.path as path
from System_tools import *


class Kraken_setup:
    rute = path.abspath(path.join(""))
    cat1 = rute + "/Cat/SCHOTT.AGF"
    cat2 = rute + "/Cat/TSPM.AGF"
    cat3 = rute + "/Cat/ESOPO.AGF"
    cat4 = rute + "/Cat/INFRARED.AGF"
    cat5 = rute + "/Cat/UTILIDADES.AGF"  # este debe de ir siempre
    filepath = [cat1, cat2, cat3, cat4, cat5]

    [CAT, NAMES, NM, ED, CD, TD, OD, LD, IT] = load_Catalog(filepath)

    file = rute + '/Cat/Alum.csv'
    W_alum, N_alum, K_alum = load_alluminum_complex(file)
