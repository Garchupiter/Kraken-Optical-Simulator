
import os.path as path
from .SystemTools import *

class Setup():
    """Kraken_setup.
    """

    rute = path.abspath(path.join(''))
    cat1 = (rute + '/KrakenOS/Cat/SCHOTT.AGF')
    cat2 = (rute + '/KrakenOS/Cat/TSPM.AGF')
    cat3 = (rute + '/KrakenOS/Cat/ESOPO.AGF')
    cat4 = (rute + '/KrakenOS/Cat/INFRARED.AGF')
    cat5 = (rute + '/KrakenOS/Cat/UTILIDADES.AGF')
    filepath = [cat1, cat2, cat3, cat4, cat5]
    [CAT, NAMES, NM, ED, CD, TD, OD, LD, IT] = load_Catalog(filepath)
    file = (rute + '/KrakenOS/Cat/Alum.csv')
    (W_alum, N_alum, K_alum) = load_alluminum_complex(file)

