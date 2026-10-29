
import numpy as np


def extract_ray_result(System):
    """Return the last traced ray data without returning the full system."""

    return {
        "nelements": System.n,
        "val": System.val,
        "Wave": System.Wave,
        "ray_SurfHits": System.ray_SurfHits,
        "SURFACE": System.SURFACE,
        "NAME": System.NAME,
        "GLASS": System.GLASS,
        "S_XYZ": System.S_XYZ,
        "T_XYZ": System.T_XYZ,
        "XYZ": System.XYZ,
        "OST_XYZ": System.OST_XYZ,
        "OST_LMN": System.OST_LMN,
        "S_LMN": System.S_LMN,
        "LMN": System.LMN,
        "R_LMN": System.R_LMN,
        "N0": System.N0,
        "N1": System.N1,
        "WAV": System.WAV,
        "G_LMN": System.G_LMN,
        "ORDER": System.ORDER,
        "GRATING": System.GRATING,
        "DISTANCE": System.DISTANCE,
        "OP": System.OP,
        "TOP_S": System.TOP_S,
        "TOP": System.TOP,
        "ALPHA": System.ALPHA,
        "BULK_TRANS": System.BULK_TRANS,
        "RP": System.RP,
        "RS": System.RS,
        "TP": System.TP,
        "TS": System.TS,
        "TTBE": System.TTBE,
        "TT": System.TT,
    }


def _as_array(seq):
    """``seq`` as an array, or as an object array when its rows differ in shape (a ray stopped part way
    keeps a placeholder record, which NumPy 2 will not stack)."""
    try:
        return np.asarray(seq)
    except ValueError:
        return np.asarray(list(seq), dtype=object)


class raykeeper():
    """raykeeper.
    """


    def __init__(self, System):
        """__init__.

        Parameters
        ----------
        System :
            System
        """
        self.SYSTEM = System
        self.clean()

    def valid(self):
        """valid.
        """
        z = np.argwhere((self.vld == 1))
        return z

    def push(self):
        """push.
        """
        self.push_result(extract_ray_result(self.SYSTEM))

    def push_result(self, result):
        """Store one traced ray from an extracted result dictionary."""
        self.nelements = result["nelements"]
        if (result["val"] == 0):
            self.invalid_vld = np.append(self.vld, 0)
            self.invalid_SURFACE.append(np.asarray(result["SURFACE"]))
            self.invalid_NAME.append(np.asarray(result["NAME"]))
            self.invalid_GLASS.append(np.asarray(result["GLASS"]))
            self.invalid_S_XYZ.append(np.asarray(result["S_XYZ"]))
            self.invalid_T_XYZ.append(np.asarray(result["T_XYZ"]))
            self.invalid_XYZ.append(np.asarray(result["XYZ"]))

            self.invalid_OST_XYZ.append(_as_array(result["OST_XYZ"]))
            self.invalid_OST_LMN.append(np.asarray(result["OST_LMN"]))
            self.invalid_S_LMN.append(np.asarray(result["S_LMN"]))
            self.invalid_LMN.append(np.asarray(result["LMN"]))
            self.invalid_R_LMN.append(np.asarray(result["R_LMN"]))
            self.invalid_N0.append(np.asarray(result["N0"]))
            self.invalid_N1.append(np.asarray(result["N1"]))
            self.invalid_WAV.append(np.asarray(result["WAV"]))
            self.invalid_G_LMN.append(np.asarray(result["G_LMN"]))
            self.invalid_ORDER.append(np.asarray(result["ORDER"]))
            self.invalid_GRATING.append(np.asarray(result["GRATING"]))
            self.invalid_DISTANCE.append(np.asarray(result["DISTANCE"]))
            self.invalid_OP.append(np.asarray(result["OP"]))
            self.invalid_TOP_S.append(np.asarray(result["TOP_S"]))
            self.invalid_TOP.append(np.asarray(result["TOP"]))
            self.invalid_ALPHA.append(_as_array(result["ALPHA"]))
            self.invalid_BULK_TRANS.append(np.asarray(result["BULK_TRANS"]))
            self.invalid_RP.append(np.asarray(result["RP"]))
            self.invalid_RS.append(np.asarray(result["RS"]))
            self.invalid_TP.append(np.asarray(result["TP"]))
            self.invalid_TS.append(np.asarray(result["TS"]))
            self.invalid_TTBE.append(np.asarray(result["TTBE"]))
            self.invalid_TT.append(np.asarray(result["TT"]))
        else:
            self.vld = np.append(self.vld, 1)
            self.valid_vld = np.append(self.vld, 0)
            self.valid_SURFACE.append(np.asarray(result["SURFACE"]))
            self.valid_NAME.append(np.asarray(result["NAME"]))
            self.valid_GLASS.append(np.asarray(result["GLASS"]))
            self.valid_S_XYZ.append(np.asarray(result["S_XYZ"]))
            self.valid_T_XYZ.append(np.asarray(result["T_XYZ"]))
            self.valid_XYZ.append(np.asarray(result["XYZ"]))
            self.valid_OST_XYZ.append(np.asarray(result["OST_XYZ"]))
            self.valid_OST_LMN.append(np.asarray(result["OST_LMN"]))
            self.valid_S_LMN.append(np.asarray(result["S_LMN"]))
            self.valid_LMN.append(np.asarray(result["LMN"]))
            self.valid_R_LMN.append(np.asarray(result["R_LMN"]))
            self.valid_N0.append(np.asarray(result["N0"]))
            self.valid_N1.append(np.asarray(result["N1"]))
            self.valid_WAV.append(np.asarray(result["WAV"]))
            self.valid_G_LMN.append(np.asarray(result["G_LMN"]))
            self.valid_ORDER.append(np.asarray(result["ORDER"]))
            self.valid_GRATING.append(np.asarray(result["GRATING"]))
            self.valid_DISTANCE.append(np.asarray(result["DISTANCE"]))
            self.valid_OP.append(np.asarray(result["OP"]))
            self.valid_TOP_S.append(np.asarray(result["TOP_S"]))
            self.valid_TOP.append(np.asarray(result["TOP"]))
            self.valid_ALPHA.append(np.asarray(result["ALPHA"]))
            self.valid_BULK_TRANS.append(np.asarray(result["BULK_TRANS"]))
            self.valid_RP.append(np.asarray(result["RP"]))
            self.valid_RS.append(np.asarray(result["RS"]))
            self.valid_TP.append(np.asarray(result["TP"]))
            self.valid_TS.append(np.asarray(result["TS"]))
            self.valid_TTBE.append(np.asarray(result["TTBE"]))
            self.valid_TT.append(np.asarray(result["TT"]))
        self.nrays = (self.nrays + 1)


        self.RayWave.append(result["Wave"])
        self.CC.append(result["ray_SurfHits"])



        self.SURFACE.append(np.asarray(result["SURFACE"]))
        self.NAME.append(np.asarray(result["NAME"]))
        self.GLASS.append(np.asarray(result["GLASS"]))
        self.S_XYZ.append(np.asarray(result["S_XYZ"]))
        self.T_XYZ.append(np.asarray(result["T_XYZ"]))
        self.XYZ.append(np.asarray(result["XYZ"]))

        # the points themselves: np.asarray(filter(None, lst)) stored a 0-d array holding the filter object
        self.OST_XYZ.append(_as_array(result["OST_XYZ"]))
        self.OST_LMN.append(np.asarray(result["OST_LMN"]))
        self.S_LMN.append(np.asarray(result["S_LMN"]))
        self.LMN.append(np.asarray(result["LMN"]))
        self.R_LMN.append(np.asarray(result["R_LMN"]))
        self.N0.append(np.asarray(result["N0"]))
        self.N1.append(np.asarray(result["N1"]))
        self.WAV.append(np.asarray(result["WAV"]))
        self.G_LMN.append(np.asarray(result["G_LMN"]))
        self.ORDER.append(np.asarray(result["ORDER"]))
        self.GRATING.append(np.asarray(result["GRATING"]))
        self.DISTANCE.append(np.asarray(result["DISTANCE"]))
        self.OP.append(np.asarray(result["OP"]))
        self.TOP_S.append(np.asarray(result["TOP_S"]))
        self.TOP.append(np.asarray(result["TOP"]))
        self.ALPHA.append(_as_array(result["ALPHA"]))
        self.BULK_TRANS.append(np.asarray(result["BULK_TRANS"]))
        self.RP.append(np.asarray(result["RP"]))
        self.RS.append(np.asarray(result["RS"]))
        self.TP.append(np.asarray(result["TP"]))
        self.TS.append(np.asarray(result["TS"]))
        self.TTBE.append(np.asarray(result["TTBE"]))
        self.TT.append(np.asarray(result["TT"]))

    def extend_results(self, results):
        """Store multiple extracted ray results."""
        for result in results:
            self.push_result(result)

    def clean(self):
        """clean.
        """
        self.vld = np.asarray([])
        self.nrays = 0
        self.RayWave = []
        self.CC =[]
        self.SURFACE = []
        self.NAME = []
        self.GLASS = []
        self.S_XYZ = []
        self.T_XYZ = []
        self.XYZ = []
        self.OST_XYZ = []
        self.OST_LMN = []
        self.S_LMN = []
        self.LMN = []
        self.R_LMN = []
        self.N0 = []
        self.N1 = []
        self.WAV = []
        self.G_LMN = []
        self.ORDER = []
        self.GRATING = []
        self.DISTANCE = []
        self.OP = []
        self.TOP_S = []
        self.TOP = []
        self.ALPHA = []
        self.BULK_TRANS = []
        self.RP = []
        self.RS = []
        self.TP = []
        self.TS = []
        self.TTBE = []
        self.TT = []
        self.valid_RayWave = []
        self.valid_CCC = []
        self.valid_SURFACE = []
        self.valid_NAME = []
        self.valid_GLASS = []
        self.valid_S_XYZ = []
        self.valid_T_XYZ = []
        self.valid_XYZ = []
        self.valid_OST_XYZ = []
        self.valid_OST_LMN = []
        self.valid_S_LMN = []
        self.valid_LMN = []
        self.valid_R_LMN = []
        self.valid_N0 = []
        self.valid_N1 = []
        self.valid_WAV = []
        self.valid_G_LMN = []
        self.valid_ORDER = []
        self.valid_GRATING = []
        self.valid_DISTANCE = []
        self.valid_OP = []
        self.valid_TOP_S = []
        self.valid_TOP = []
        self.valid_ALPHA = []
        self.valid_BULK_TRANS = []
        self.valid_RP = []
        self.valid_RS = []
        self.valid_TP = []
        self.valid_TS = []
        self.valid_TTBE = []
        self.valid_TT = []
        self.invalid_RayWave = []
        self.invalid_CCC = []
        self.invalid_SURFACE = []
        self.invalid_NAME = []
        self.invalid_GLASS = []
        self.invalid_S_XYZ = []
        self.invalid_T_XYZ = []
        self.invalid_XYZ = []
        self.invalid_OST_XYZ = []
        self.invalid_OST_LMN = []
        self.invalid_S_LMN = []
        self.invalid_LMN = []
        self.invalid_R_LMN = []
        self.invalid_N0 = []
        self.invalid_N1 = []
        self.invalid_WAV = []
        self.invalid_G_LMN = []
        self.invalid_ORDER = []
        self.invalid_GRATING = []
        self.invalid_DISTANCE = []
        self.invalid_OP = []
        self.invalid_TOP_S = []
        self.invalid_TOP = []
        self.invalid_ALPHA = []
        self.invalid_BULK_TRANS = []
        self.invalid_RP = []
        self.invalid_RS = []
        self.invalid_TP = []
        self.invalid_TS = []
        self.invalid_TTBE = []
        self.invalid_TT = []

    def pick(self, N_ELEMENT=(- 1), coordinates = "global"):
        """pick.

        Parameters
        ----------
        N_ELEMENT :
            N_ELEMENT

            coordinates = "global" or "local"
        """

        gls = self.SYSTEM.SDT[N_ELEMENT].Glass
        if gls == "NULL":
            print("NULL surface has been chosen, the return values correspond to those of the previous surface")

        self.numsup = (self.nelements - 1)

        if coordinates == "global":
            self.xyz = self.valid_XYZ
            self.lmn = self.valid_LMN
        else:
            self.xyz = self.valid_OST_XYZ
            self.lmn = self.valid_OST_LMN

        self.s = self.valid_SURFACE
        if ((N_ELEMENT < 0) or (N_ELEMENT > self.numsup)):
            N_ELEMENT = self.numsup
        else:
            N_ELEMENT = N_ELEMENT
        # AA = []
        BB = []
        for k in self.s:
            aa = np.argwhere((k == N_ELEMENT))
            aa = np.squeeze(aa)
            # print(aa)
            # AA.append(aa)
            BB.append(np.size(aa))
        # AA = np.asarray(AA)
        BB = np.asarray(BB)
        if (N_ELEMENT != 0):
            BB = np.argwhere((BB == 1))
        else:
            BB = np.argwhere((BB == 0))
        X = []
        Y = []
        Z = []
        L = []
        M = []
        N = []
        for c in BB:
            for d in c:
                ray0 = self.xyz[d]

                [x1, y1, z1] = ray0[N_ELEMENT]
                X.append(x1)
                Y.append(y1)
                Z.append(z1)
                ray1 = self.lmn[d]
                if (N_ELEMENT != 0):
                    el = (N_ELEMENT - 1)
                else:
                    el = 0
                [l1, m1, n1] = ray1[el]
                L.append(l1)
                M.append(m1)
                N.append(n1)
        return (np.asarray(X), np.asarray(Y), np.asarray(Z), np.asarray(L), np.asarray(M), np.asarray(N))
