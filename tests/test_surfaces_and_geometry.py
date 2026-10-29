"""Surface, geometry and tracing regressions.

- Each surface has its own default arrays and lists (they were made once, in the signature, and
  shared: s.ZNK[3] = 0.1 on one surface changed every surface built without its own).
- RvTrace from the stop retraces the forward ray (it started two media off), and NsTrace's NULL
  record has the fields __CollectData unpacks (it left out LMNObjSpace).
- SetData and SetSolid keep the IgnoreVignetting setting (they made a fresh InterNormalCalc).
- cat2surf gives 200 aspheric terms (a shorter array failed in Kos.system), the block's last
  thickness, the tilts with AxisMove 0, and a first surface with any other AxisMove.
"""

import contextlib
import inspect
import io
import re

import numpy as np


def _system(surfaces, build=1):
    import KrakenOS as Kos

    with contextlib.redirect_stdout(io.StringIO()):
        return Kos.system(surfaces, Kos.Setup(), build=build)


def _surf(**values):
    import KrakenOS as Kos

    surface = Kos.surf()
    for name, value in values.items():
        setattr(surface, name, value)
    return surface


def test_surfaces_do_not_share_their_default_arrays():
    import KrakenOS as Kos

    first, second = Kos.surf(), Kos.surf()
    first.ZNK[3] = 0.1
    first.AspherData[1] = 1e-6
    first.Coating[0].append(1.0)
    first.Const[1] = 1
    third = Kos.surf()
    for other in (second, third):
        assert other.ZNK[3] == 0.0 and other.AspherData[1] == 0.0
        assert other.Coating[0] == [] and other.Const[1] == 0


def _reverse_system():
    return _system([
        _surf(Thickness=50.0, Diameter=40.0),
        _surf(Rc=40.0, Thickness=5.0, Glass="BK7", Diameter=30.0),
        _surf(Rc=-40.0, Thickness=10.0, Diameter=30.0),
        _surf(Thickness=10.0, Diameter=20.0),                      # the stop
        _surf(Rc=30.0, Thickness=4.0, Glass="F2", Diameter=30.0),
        _surf(Rc=-60.0, Thickness=40.0, Diameter=30.0),
        _surf(Diameter=80.0),
    ])


def test_a_ray_traced_back_from_the_stop_retraces_the_forward_ray():
    system = _reverse_system()
    direction = np.array([0.0, 0.05, 1.0])
    direction /= np.linalg.norm(direction)
    with contextlib.redirect_stdout(io.StringIO()):
        system.Trace([0.0, 2.0, 0.0], list(direction), 0.55)
        forward = [np.asarray(p, dtype=float) for p in system.XYZ[:4]]      # object, 1, 2, stop
        arriving = np.asarray(system.LMN[2], dtype=float)
        system.RvTrace(list(forward[3]), list(-arriving), 0.55, 3)
    back = [np.asarray(p, dtype=float) for p in system.XYZ]
    assert [int(v) for v in system.SURFACE] == [3, 2, 1, 0]
    for k in range(4):
        assert np.allclose(back[k + 1], forward[3 - k], atol=1e-6), (k, back[k + 1], forward[3 - k])


def test_the_null_surface_record_has_the_fields_collect_data_unpacks():
    from KrakenOS import KrakenSys

    source = inspect.getsource(KrakenSys.system.NsTrace)
    null_record = re.search(r"if self\.Glass\[j\] == 'NULL':\s*\n\s*ang = 0\s*\n(?:\s*#.*\n)*\s*ValToSav = \[([^\]]*)\]", source)
    unpacked = re.search(r"\[([^\]]*)\] = ValToSav", inspect.getsource(KrakenSys.system._system__CollectData))
    fields = lambda text: [name.strip() for name in text.split(",")]
    assert null_record and unpacked
    assert fields(null_record.group(1)) == fields(unpacked.group(1))


def test_set_data_keeps_the_vignetting_setting():
    system = _reverse_system()
    system.IgnoreVignetting()
    with contextlib.redirect_stdout(io.StringIO()):
        system.SetData()
    assert (system.INORM.Disable_Inner, system.INORM.ExtraDiameter) == (0, 1)


LENS = {
    "SUFR 1": {"Rc": 1 / 40.0, "Thickness": 0.0, "Glass": "BK7", "Diameter": 25.0, "aspherics": [0.0, 2e-6]},
    "SUFR 2": {"Rc": -1 / 40.0, "Thickness": 5.0, "Glass": "AIR", "Diameter": 25.0},
}


def test_a_catalog_lens_becomes_surfaces_a_system_builds():
    from KrakenOS.LensCat import cat2surf

    surfaces = cat2surf(LENS, Thickness=12.0)
    assert [len(np.asarray(s.AspherData)) for s in surfaces] == [200, 200]
    assert float(surfaces[0].AspherData[1]) == 2e-6
    _system([_surf(Thickness=50.0, Diameter=25.0), *surfaces, _surf(Diameter=60.0)])
    assert float(surfaces[-1].Thickness) == 12.0
    assert float(cat2surf(LENS, Thickness=12.0, inverse=True)[-1].Thickness) == 12.0


def test_a_catalog_lens_is_placed_as_asked():
    from KrakenOS.LensCat import cat2surf

    assert [float(s.TiltX) for s in cat2surf(LENS, TiltX=2.0, AxisMove=0)] == [2.0, 2.0]
    twice = cat2surf(LENS, TiltX=2.0, AxisMove=2)
    assert float(twice[0].TiltX) == 2.0 and float(twice[-1].AxisMove) == 2.0
