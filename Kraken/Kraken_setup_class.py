# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 17:55:21 2021

@author: JOELHERRERAVAZQUEZ
"""

import os.path as path
from .System_tools import *


class Kraken_setup:
    rute = path.abspath(path.join(""))
    cat1 = rute + "/Kraken/Cat/SCHOTT.AGF"
    cat2 = rute + "/Kraken/Cat/TSPM.AGF"
    cat3 = rute + "/Kraken/Cat/ESOPO.AGF"
    cat4 = rute + "/Kraken/Cat/INFRARED.AGF"
    cat5 = rute + "/Kraken/Cat/UTILIDADES.AGF"  # este debe de ir siempre
    filepath = [cat1, cat2, cat3, cat4, cat5]

    [CAT, NAMES, NM, ED, CD, TD, OD, LD, IT] = load_Catalog(filepath)

    file = rute + '/Kraken/Cat/Alum.csv'
    W_alum, N_alum, K_alum = load_alluminum_complex(file)
