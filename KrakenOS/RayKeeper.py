
import numpy as np


_ARRAY_FIELDS = (
    "SURFACE",
    "NAME",
    "GLASS",
    "S_XYZ",
    "T_XYZ",
    "XYZ",
    "OST_XYZ",
    "OST_LMN",
    "S_LMN",
    "LMN",
    "R_LMN",
    "N0",
    "N1",
    "WAV",
    "G_LMN",
    "ORDER",
    "GRATING",
    "DISTANCE",
    "OP",
    "TOP_S",
    "TOP",
    "ALPHA",
    "BULK_TRANS",
    "RP",
    "RS",
    "TP",
    "TS",
    "TTBE",
    "TT",
)

_FILTERED_FIELDS = frozenset(("OST_XYZ", "ALPHA"))


def _copy_field(value):
    if isinstance(value, np.ndarray):
        return value.copy()
    if isinstance(value, list):
        return [_copy_field(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_copy_field(item) for item in value)
    return value


def _is_filter_none_keep(item):
    """Match historical filter(None, ...) without NumPy truth-value errors."""
    if item is None or item is False:
        return False
    if isinstance(item, np.ndarray):
        return item.size > 0
    try:
        return bool(item)
    except (ValueError, TypeError):
        return True


def _as_stored_array(value, filtered=False):
    if filtered:
        try:
            value = [item for item in value if _is_filter_none_keep(item)]
        except TypeError:
            pass
    return np.asarray(value)


def extract_ray_result(System, copy=False):
    """Return the last traced ray data without returning the full system.

    Parameters
    ----------
    System :
        Optical system after Trace/NsTrace.
    copy : bool, optional
        If True, return a snapshot safe for retention or IPC.
        If False (default), return shallow references for immediate push().
    """

    result = {
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
    if copy:
        return {
            key: _copy_field(value) if key not in ("nelements", "val", "Wave") else value
            for key, value in result.items()
        }
    return result


class _BucketView:
    """Read-only view over the universal ray lists filtered by validity."""

    __slots__ = ("_owner", "_attr", "_want_valid")

    def __init__(self, owner, attr, want_valid):
        self._owner = owner
        self._attr = attr
        self._want_valid = want_valid

    def _indices(self):
        want = self._want_valid
        return [index for index, flag in enumerate(self._owner._ray_valid) if flag is want]

    def __len__(self):
        want = self._want_valid
        return sum(1 for flag in self._owner._ray_valid if flag is want)

    def __getitem__(self, key):
        items = getattr(self._owner, self._attr)
        indices = self._indices()
        if isinstance(key, slice):
            return [items[index] for index in indices[key]]
        return items[indices[key]]

    def __iter__(self):
        items = getattr(self._owner, self._attr)
        want = self._want_valid
        for flag, item in zip(self._owner._ray_valid, items):
            if flag is want:
                yield item

    def __repr__(self):
        return f"<_BucketView attr={self._attr!r} want_valid={self._want_valid} len={len(self)}>"


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

    @property
    def vld(self):
        if self._vld_cache is None:
            self._vld_cache = np.asarray(self._vld_list, dtype=float)
        return self._vld_cache

    def _invalidate_vld_cache(self):
        self._vld_cache = None

    def valid(self):
        """valid.
        """
        z = np.argwhere((self.vld == 1))
        return z

    def push(self):
        """push.
        """
        self.push_result(extract_ray_result(self.SYSTEM, copy=False))

    def push_result(self, result):
        """Store one traced ray from an extracted result dictionary."""
        self.nelements = result["nelements"]
        is_valid = bool(result["val"] != 0)

        stored = {}
        for field in _ARRAY_FIELDS:
            stored[field] = _as_stored_array(
                result[field],
                filtered=(field in _FILTERED_FIELDS),
            )

        if is_valid:
            self._vld_list.append(1)
            self._invalidate_vld_cache()
            self.valid_vld = np.asarray(self._vld_list + [0], dtype=float)
        else:
            self.invalid_vld = np.asarray(self._vld_list + [0], dtype=float)

        self._ray_valid.append(is_valid)
        self.nrays = self.nrays + 1

        self.RayWave.append(result["Wave"])
        hits = result["ray_SurfHits"]
        if isinstance(hits, np.ndarray):
            self.CC.append(hits)
        else:
            self.CC.append(np.asarray(hits))

        for field in _ARRAY_FIELDS:
            getattr(self, field).append(stored[field])

    def extend_results(self, results):
        """Store multiple extracted ray results."""
        for result in results:
            self.push_result(result)

    def _make_bucket_views(self):
        for field in _ARRAY_FIELDS:
            setattr(self, "valid_" + field, _BucketView(self, field, True))
            setattr(self, "invalid_" + field, _BucketView(self, field, False))

    def clean(self):
        """clean.
        """
        self._vld_list = []
        self._vld_cache = None
        self._ray_valid = []
        self.valid_vld = np.asarray([], dtype=float)
        self.invalid_vld = np.asarray([], dtype=float)
        self.nelements = 0
        self.nrays = 0
        self.RayWave = []
        self.CC = []
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
        self.invalid_RayWave = []
        self.invalid_CCC = []
        self._make_bucket_views()
        for attr in ("numsup", "xyz", "lmn", "s"):
            if hasattr(self, attr):
                delattr(self, attr)

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
