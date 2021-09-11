# PIP
import site
import os.path as path
from .SystemTools import *

import os.path # comentar para local
class Setup():
    """Kraken_setup.
    """

    rute0 = site.getsitepackages()
    rute = rute0[0]

    print("-----------------------------")

    print("  The library is installed in")
    print(rute + '/KrakenOS/')
    print("  For help see Examples directory")
    print("  To change glass catalog modify catalog list in SetupCas.py")
    print("  See User_Manual USER_MANUAL_KrakenOS_Provisional.pdf in:" )
    print(print(rute + '/Docs/'))
    print("-----------------------------")

    ## Add catalog you need
    cat1 = (rute + '/KrakenOS/Cat/SCHOTT.AGF')
    cat2 = (rute + '/KrakenOS/Cat/TSPM.AGF')
    cat3 = (rute + '/KrakenOS/Cat/ESOPO.AGF')
    cat4 = (rute + '/KrakenOS/Cat/INFRARED.AGF')
    cat5 = (rute + '/KrakenOS/Cat/UTILIDADES.AGF')

    ## Add catalog to the list filepath
    filepath = [cat1, cat2, cat3, cat4, cat5]

    [CAT, NAMES, NM, ED, CD, TD, OD, LD, IT] = load_Catalog(filepath)
    file = (rute + '/KrakenOS/Cat/Alum.csv')

    (W_alum, N_alum, K_alum) = load_alluminum_complex(file)

