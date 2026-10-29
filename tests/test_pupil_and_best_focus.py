"""Pupil and best-focus regressions.

- BestFocus returns the shift of least RMS radius; it asked fsolve for a zero of the RMS, which a
  real bundle has not, and returned wherever the solver stopped.
- Fpupil reads SYSTEM.SURFACE (it read SurfACE) and the hit on the surface it is given.
- RayStop sets the launch direction its calculateXYZ reads (it set IMP_LMN, read INP_LMN).
- PupilCalc finds a stop with no gap before it (an object at infinity with the stop as the first
  surface): its Newton step was 0/0 at r = 0 and it raised IndexError. The rays it then launches
  leave along the field instead of 0/0. A telecentric object space says why it has no pupil.
- Phase's off-axis wavefront with such a stop is the one the same lens gives with a gap before the
  stop (an object at infinity: the gap changes nothing). Phase's input plane, tilted through the
  chief ray's start, crossed the stop, and half the pupil's paths were counted the wrong way.
"""

import contextlib
import io

import numpy as np
import pytest


def _surf(Kos, **values):
    surface = Kos.surf()
    for name, value in values.items():
        setattr(surface, name, value)
    return surface


def _system(surfaces):
    import KrakenOS as Kos

    with contextlib.redirect_stdout(io.StringIO()):
        return Kos.system(surfaces, Kos.Setup(), build=1)


def _singlet(gap_before_lens=10.0):
    import KrakenOS as Kos

    return _system([
        _surf(Kos, Thickness=gap_before_lens, Diameter=30.0),
        _surf(Kos, Rc=50.0, Thickness=6.0, Glass="BK7", Diameter=30.0),
        _surf(Kos, Rc=-50.0, Thickness=45.0, Diameter=30.0),
        _surf(Kos, Diameter=80.0),
    ])


def _r_rms(delta, L, M, N, X, Y):
    from KrakenOS.LibRMS import R_RMS_delta

    return R_RMS_delta(delta, L, M, N, X, Y)


def _scan(L, M, N, X, Y, around):
    grid = np.linspace(around - 30.0, around + 30.0, 60001)
    values = [_r_rms(dz, L, M, N, X, Y) for dz in grid]
    return float(grid[int(np.argmin(values))])


def test_best_focus_of_an_astigmatic_bundle_is_its_least_rms():
    from KrakenOS.LibRMS import BestFocus

    rng = np.random.default_rng(1111)
    u, v = rng.uniform(-0.05, 0.05, 400), rng.uniform(-0.05, 0.05, 400)
    N = 1.0 / np.sqrt(1.0 + u * u + v * v)
    L, M = u * N, v * N
    X, Y = 5.0 * u, 15.0 * v                       # x focuses at dz = -5, y at dz = -15
    best = BestFocus(X, Y, np.zeros_like(X), L, M, N)
    assert best == pytest.approx(_scan(L, M, N, X, Y, -10.0), abs=2e-3)
    assert _r_rms(best, L, M, N, X, Y) < _r_rms(best - 0.1, L, M, N, X, Y)
    assert _r_rms(best, L, M, N, X, Y) < _r_rms(best + 0.1, L, M, N, X, Y)


def test_best_focus_of_a_singlet_bundle():
    import KrakenOS as Kos
    from KrakenOS.LibRMS import BestFocus

    system = _singlet(1000.0)
    rays = Kos.raykeeper(system)
    for radius in np.linspace(0.5, 4.0, 8):
        for angle in np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False):
            system.Trace([radius * np.cos(angle), radius * np.sin(angle), 0.0], [0.0, 0.0, 1.0], 0.55)
            rays.push()
    x, y, z, l, m, n = rays.pick(-1)
    best = BestFocus(x, y, z, l, m, n)
    assert abs(best) > 0.1
    assert best == pytest.approx(_scan(l, m, n, x, y, 0.0), abs=1e-3)


def test_best_focus_of_parallel_rays_is_zero():
    from KrakenOS.LibRMS import BestFocus

    assert BestFocus(np.linspace(-1, 1, 9), np.zeros(9), np.zeros(9), np.zeros(9), np.zeros(9), np.ones(9)) == 0.0


def test_fpupil_measures_the_ray_where_it_meets_the_surface():
    from KrakenOS.PupilTool import Fpupil

    residual = float(Fpupil(_singlet(), 3.0, 1.0, 1, 0.0, 0.55, 1))
    assert residual == pytest.approx(2.0, abs=1e-6)      # a ray at y = 3 mm, target 1 mm


def test_raystop_has_the_launch_direction_it_traces_with():
    from KrakenOS.PupilTool import RayStop

    with contextlib.redirect_stdout(io.StringIO()):
        stop = RayStop(_singlet(), 0.55, 1)
        stop.Desp_X = stop.Desp_Y = 0.0
        landed = [float(v) for v in stop.calculateXYZ(0.0, 2.0)]
    assert list(stop.INP_LMN) == [0, 0, 1]
    assert np.all(np.isfinite(landed))


def _gapless_stop_system(gap=0.0):
    """An object at infinity, the stop right after it, then a 50 mm lens (the F-theta layout's case)."""
    import KrakenOS as Kos

    return _system([
        _surf(Kos, Thickness=gap, Diameter=10.0),
        _surf(Kos, Thickness=5.0, Diameter=2.0),
        _surf(Kos, Rc=51.5, Thickness=4.0, Glass="BK7", Diameter=20.0),
        _surf(Kos, Rc=-51.5, Thickness=48.4, Diameter=20.0),
        _surf(Kos, Diameter=60.0),
    ])


def _pupil(system, stop):
    import KrakenOS as Kos

    with contextlib.redirect_stdout(io.StringIO()):
        return Kos.PupilCalc(system, stop, 0.55, "EPD", 2.0)


def test_pupil_of_a_stop_with_no_gap_before_it():
    pupil = _pupil(_gapless_stop_system(), 1)
    assert float(np.ravel(pupil.PosPupInp)[2]) == pytest.approx(0.0, abs=1e-9)
    assert float(np.ravel(pupil.RadPupInp)[0]) == pytest.approx(1.0)
    pupil.Samp, pupil.Ptype, pupil.FieldType, pupil.FieldX = 6, "hexapolar", "angle", 0.0
    for field in (0.0, 10.0):
        pupil.FieldY = field
        _x, _y, _z, L, M, N = pupil.Pattern2Field()
        assert np.all(np.isfinite(L)) and np.all(np.isfinite(M)) and np.all(np.isfinite(N))
        assert np.allclose(M, np.sin(np.deg2rad(field))) and np.allclose(N, np.cos(np.deg2rad(field)))


def test_pupil_of_a_stop_after_a_gap_is_unchanged():
    import KrakenOS as Kos

    system = _system([
        _surf(Kos, Thickness=10.0, Diameter=10.0),
        _surf(Kos, Thickness=5.0, Diameter=2.0),
        _surf(Kos, Rc=51.5, Thickness=4.0, Glass="BK7", Diameter=20.0),
        _surf(Kos, Rc=-51.5, Thickness=48.4, Diameter=20.0),
        _surf(Kos, Diameter=60.0),
    ])
    assert float(np.ravel(_pupil(system, 1).PosPupInp)[2]) == pytest.approx(10.0, abs=1e-6)


def test_a_telecentric_object_space_has_no_pupil_to_find():
    import KrakenOS as Kos

    system = _system([
        _surf(Kos, Thickness=100.0, Diameter=30.0),
        _surf(Kos, Thickness=50.0, Diameter=30.0, Thin_Lens=50.0),
        _surf(Kos, Thickness=100.0, Diameter=5.0),            # the stop, at the lens's back focal point
        _surf(Kos, Diameter=60.0),
    ])
    with pytest.raises(ValueError, match="infinity"):
        _pupil(system, 2)


def _wavefront(system, field):
    import KrakenOS as Kos

    pupil = _pupil(system, 1)
    pupil.Samp, pupil.Ptype, pupil.FieldType, pupil.FieldX, pupil.FieldY = 10, "hexapolar", "angle", 0.0, field
    with contextlib.redirect_stdout(io.StringIO()):
        y, x, w, _p2v = Kos.Phase(pupil)
    return np.asarray(x), np.asarray(y), np.asarray(w)


def test_phase_off_axis_with_no_gap_before_the_stop():
    for field in (5.0, 10.0):
        _x, _y, with_gap = _wavefront(_gapless_stop_system(10.0), field)
        _x, _y, without = _wavefront(_gapless_stop_system(0.0), field)
        assert np.ptp(with_gap) < 5.0, np.ptp(with_gap)        # the bug gave hundreds of waves
        assert np.ptp(without) == pytest.approx(np.ptp(with_gap), rel=1e-3, abs=1e-4)


def test_the_exit_pupil_searches_find_the_least_rms_too():
    from KrakenOS.LibRMS import BestFocus
    from KrakenOS.PhaseCalc import BestExitPupilPos as phase_exit_pupil
    from KrakenOS.SphericalPSFCalc import BestExitPupilPos as psf_exit_pupil

    rng = np.random.default_rng(1119)
    u, v = rng.uniform(-0.05, 0.05, 400), rng.uniform(-0.05, 0.05, 400)
    N = 1.0 / np.sqrt(1.0 + u * u + v * v)
    L, M = u * N, v * N
    X, Y, Z = 5.0 * u, 15.0 * v, np.zeros(400)
    best = BestFocus(X, Y, Z, L, M, N)
    assert best == pytest.approx(_scan(L, M, N, X, Y, -10.0), abs=2e-3)
    assert phase_exit_pupil(X, Y, Z, L, M, N) == pytest.approx(best)
    assert psf_exit_pupil(X, Y, Z, L, M, N) == pytest.approx(best)
