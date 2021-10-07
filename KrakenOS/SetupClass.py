""" change pip to True for pip installation """


import inspect
from .SystemTools import *
RUTE=inspect.getmodule(load_Catalog).__file__
rute=RUTE[:-15]+ "/Cat/"





class Setup():

    def __init__(self):
        """Kraken_setup.
        """
        print("Default catalog is loaded from: ", rute)
        cat1 = (rute + 'SCHOTT.AGF')
        cat2 = (rute + 'TSPM.AGF')
        cat3 = (rute + 'ESOPO.AGF')
        cat4 = (rute + 'INFRARED.AGF')

        self.GlassCat =[]
        self.GlassCat.append(cat1)
        self.GlassCat.append(cat2)
        self.GlassCat.append(cat3)
        self.GlassCat.append(cat4)

        self.Load(self.GlassCat)

    def Load(self, GL):


        utilities = (rute + 'UTILIDADES.AGF')
        GL.append(utilities)


        [self.CAT, self.NAMES, self.NM, self.ED, self.CD, self.TD, self.OD, self.LD, self.IT] = load_Catalog(GL)
        file = (rute + 'Alum.csv')
        (self.W_alum, self.N_alum, self.K_alum) = load_alluminum_complex(file)




