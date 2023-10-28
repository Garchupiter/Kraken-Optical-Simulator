import pandas as pd
import inspect
import os
from .SystemTools import *
RUTE=inspect.getmodule(load_Catalog).__file__
rute=RUTE[:-15]+ "/Cat/"


def read_material_data(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        section = None
        data = {
            'wl': [],
            'n': [],
            'k': [],
        }

        for line in lines:
            line = line.strip()
            if line == 'wl,n':
                section = 'n'
            elif line == 'wl,k':
                section = 'k'
            elif line:
                if section == 'n':
                    wl, n = map(float, line.split(','))
                    data['wl'].append(wl)
                    data['n'].append(n)
                elif section == 'k':
                    wl, k = map(float, line.split(','))
                    # data['wl'].append(wl)  # Se añade la longitud de onda a ambos n y k
                    data['k'].append(k)

        w = np.array(data['wl'])
        n = np.array(data['n'])
        k = np.array(data['k'])

        # print(data['wl'])

        return w,n,k

    except Exception as e:
        print(f"Error: {e}")
        return None





class Setup():

    def __init__(self):
        """Kraken_setup.
        """
        print("Default catalog is loaded from: ", rute)

        # cat1 = (rute + 'SCHOTT.AGF')
        # cat2 = (rute + 'infrared.agf')

        self.GlassCat =[rute + cat for cat in os.listdir(rute) if cat.endswith(('.AGF', '.agf'))]
        # self.GlassCat.append(cat1)
        # self.GlassCat.append(cat2)

        self.Load(self.GlassCat)

        file = (rute + 'Alum.csv')
        self.W_met, self.N_met, self.K_met, self.Name_met = [], [], [], []
        (W_alum, N_alum, K_alum) = load_metal_complex(file)
        self.W_met.append(W_alum)
        self.N_met.append(N_alum)
        self.K_met.append(K_alum)
        self.Name_met.append("Alum")

    def Load(self, GL):

        print("Glass catalog loaded")
        utilities = (rute + 'UTILIDADES.AGF')
        GL.append(utilities)

        [self.CAT, self.NAMES, self.NM, self.ED, self.CD, self.TD, self.OD, self.LD, self.IT] = load_Catalog(GL)


    def LoadMetal(self, FILE, nm = "Name", Type = 1):
        if Type == 0:
            (W_met, N_met, K_met) = load_metal_complex(FILE)
            self.W_met.append(W_met)
            self.N_met.append(N_met)
            self.K_met.append(K_met)
            self.Name_met.append(nm)

        if Type == 1:

            # Lee el archivo CSV
            Wl, n, k = read_material_data(FILE)

            self.W_met.append(Wl)
            self.N_met.append(n)
            self.K_met.append(k)
            self.Name_met.append(nm)
