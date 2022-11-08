import inspect
from .SystemTools import *
RUTE=inspect.getmodule(load_Catalog).__file__
rute=RUTE[:-15]+ "/Cat/"

class Setup():

    def __init__(self):
        """Kraken_setup.
        """
        print("Catalog dir: ", rute)
        cat1 = (rute + 'SCHOTT.AGF')
        cat2 = (rute + 'infrared.agf')

        self.GlassCat =[]
        self.GlassCat.append(cat1)
        self.GlassCat.append(cat2)

        self.Load(self.GlassCat)

    def Load(self, GL):

        print("Glass catalog loaded")
        utilities = (rute + 'UTILIDADES.AGF')
        GL.append(utilities)

        [self.CAT, self.NAMES, self.NM, self.ED, self.CD, self.TD, self.OD, self.LD, self.IT] = load_Catalog(GL)
        file = (rute + 'Alum.csv')
        (self.W_alum, self.N_alum, self.K_alum) = load_alluminum_complex(file)




