"""raykeeper stores OST_XYZ and ALPHA as the trace's values for every ray.

push stored them as np.asarray(filter(None, lst)) -- in Python 3 a 0-d object array holding the
filter object, not the points -- in rays.OST_XYZ / rays.ALPHA and, for an invalid ray, in
rays.invalid_OST_XYZ / rays.invalid_ALPHA. (valid_OST_XYZ / valid_ALPHA already held the lists.)
"""

import numpy as np

from test_raykeeper_results import build_simple_system, trace_invalid_ray, trace_valid_ray


def _numeric(value):
    array = np.asarray(value)
    return array.ndim >= 1 and array.dtype != object and np.all(np.isfinite(array.astype(float)))


def test_a_valid_ray_keeps_its_object_space_points_and_alphas():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    trace_valid_ray(system)
    rays.push()
    assert _numeric(rays.OST_XYZ[0]) and np.asarray(rays.OST_XYZ[0]).shape[-1] == 3
    assert _numeric(rays.ALPHA[0])
    assert np.allclose(np.asarray(rays.OST_XYZ[0], dtype=float), np.asarray(rays.valid_OST_XYZ[0], dtype=float))


def test_an_invalid_ray_keeps_its_object_space_points_and_alphas():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    trace_invalid_ray(system)
    rays.push()
    for stored in (rays.OST_XYZ[0], rays.ALPHA[0], rays.invalid_OST_XYZ[0], rays.invalid_ALPHA[0]):
        assert not isinstance(np.asarray(stored).reshape(-1)[0] if np.asarray(stored).size else None, filter)
        assert np.asarray(stored).ndim >= 1
