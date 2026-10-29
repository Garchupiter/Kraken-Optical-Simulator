"""Wavefront, PSF and MTF regressions.

- Phase leaves the caller's system as it was: its surfaces keep their masks (EraseVTK ran on them
  before they were copied, and the copies kept Mask_Type with no shape, so a masked system's
  wavefront raised) and the system its vignetting setting.
- Zernike_Fitting runs under NumPy 2 (it set 1 x 1 matrices into a float array), and with an
  integer selector (np.ones(n, dtype=int)) returns the coefficients, not integers.
- psf returns its extent when asked; PSFMap.generate_psf_map unpacked one it never returned.
- calculate_mtf uses its pixels and PupilSample (it always used 600 and 4); psf4mtf and psf take a
  perfect lens's all-zero coefficients; plot_mtf, given the focal length, puts a perfect lens's
  MTF cutoff at Diameter / (wavelength * Focal).
"""

import contextlib
import io
import os

import numpy as np
import pytest


def _masked_singlet():
    import pyvista as pv

    import KrakenOS as Kos

    obj = Kos.surf(); obj.Thickness = 1000.0; obj.Diameter = 30.0
    front = Kos.surf(); front.Rc = 50.0; front.Thickness = 6.0; front.Glass = "BK7"; front.Diameter = 30.0
    shape = pv.MultiBlock()
    shape.append(pv.Disc(center=[20.0, 0.0, 0.0], inner=0, outer=1.0, normal=[0, 0, 1], r_res=1, c_res=60))
    front.Mask_Shape = shape
    front.Mask_Type = 2                       # an obstruction outside the beam: nothing is stopped
    back = Kos.surf(); back.Rc = -50.0; back.Thickness = 47.0; back.Diameter = 30.0
    image = Kos.surf(); image.Diameter = 60.0
    with contextlib.redirect_stdout(io.StringIO()):
        return Kos.system([obj, front, back, image], Kos.Setup(), build=1)


def test_phase_leaves_the_callers_system_as_it_was():
    import KrakenOS as Kos

    system = _masked_singlet()
    before = (system.INORM.Disable_Inner, system.INORM.ExtraDiameter)
    with contextlib.redirect_stdout(io.StringIO()):
        pupil = Kos.PupilCalc(system, 1, 0.55, "EPD", 8.0)
        pupil.Samp, pupil.Ptype, pupil.FieldType, pupil.FieldX, pupil.FieldY = 6, "hexapolar", "angle", 0.0, 0.0
        _x, _y, wave, p2v = Kos.Phase(pupil)
    assert not isinstance(system.SDT[1].Mask_Shape, str)            # EraseVTK sets it to "None"
    assert (system.INORM.Disable_Inner, system.INORM.ExtraDiameter) == before
    wave = np.asarray(wave, dtype=float)
    assert wave.size > 10 and np.all(np.isfinite(wave)) and 0.0 < float(p2v) < 5.0


def test_zernike_fit_with_an_integer_selector():
    from KrakenOS.WavefrontFit import Zernike_Fitting

    rng = np.random.default_rng(1120)
    x, y = rng.uniform(-0.7, 0.7, 400), rng.uniform(-0.7, 0.7, 400)
    z = 0.3 * (2.0 * (x * x + y * y) - 1.0) + 0.1 * x
    as_floats = np.asarray(Zernike_Fitting(x, y, z, np.ones(15))[0], dtype=float)
    as_ints = np.asarray(Zernike_Fitting(x, y, z, np.ones(15, dtype=int))[0], dtype=float)
    assert np.any(np.abs(as_floats) > 1e-3)
    assert np.allclose(as_ints, as_floats)


def test_psf_returns_its_extent_when_asked():
    from KrakenOS.PSFCalc import psf

    coef = np.zeros(15)
    coef[4] = 0.1
    alone = psf(coef, 50.0, 8.0, 0.55, pixels=64, PupilSample=4)
    image, extent = psf(coef, 50.0, 8.0, 0.55, pixels=64, PupilSample=4, return_extent=True)
    assert np.array_equal(alone, image)
    assert len(extent) == 4 and extent[1] > extent[0]


def test_psf_map_makes_its_map(tmp_path):
    import KrakenOS as Kos

    obj = Kos.surf(); obj.Thickness = 1000.0; obj.Diameter = 30.0
    front = Kos.surf(); front.Rc = 50.0; front.Thickness = 6.0; front.Glass = "BK7"; front.Diameter = 30.0
    back = Kos.surf(); back.Rc = -50.0; back.Thickness = 47.0; back.Diameter = 30.0
    image = Kos.surf(); image.Diameter = 60.0
    cwd = os.getcwd()
    os.chdir(tmp_path)                       # generate_psf_map saves numpyPSFs<wave>.npy where it runs
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            system = Kos.system([obj, front, back, image], Kos.Setup(), build=1)
            mapper = Kos.PSFMap(system, fieldxs=[0, 0], fieldys=[0, 1], surf=1, aperVal=8.0, aperType="EPD",
                                stepx=1, stepy=1, Wx=0.55)
            mapper.calc_wide_angle_pupil_list(sampling=3)
            # the wide-angle path: calc_wide_angle_pupil_list makes single-ray ("rtheta") pupils, which
            # the default wideAngle=False path hands to Phase and cannot fit
            mapper.generate_psf_map(pixels=32, PupilSample=4, plot=0, NC=15, wideAngle=True)
    finally:
        os.chdir(cwd)
    assert len(mapper.PSFList) >= 1


def test_mtf_of_a_perfect_lens_ends_at_its_cutoff():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from KrakenOS.PSFCalc import calculate_mtf, plot_mtf

    focal, diameter, wave, sample = 50.0, 8.0, 0.55, 4
    mtf = np.asarray(calculate_mtf(np.zeros(15), focal, diameter, wave, pixels=256, PupilSample=sample))
    assert mtf.shape == (256, 256)
    plt.close("all")
    plot_mtf(mtf, diameter, wave, Focal=focal, PupilSample=sample)
    line = plt.gca().get_lines()[0]
    freq, value = np.asarray(line.get_xdata()), np.asarray(line.get_ydata())
    plt.close("all")
    cutoff = diameter / (wave * 1e-3 * focal)
    step = freq[1] - freq[0]
    assert abs(freq[np.nonzero(value > 1e-6)[0][-1]] - cutoff) <= 1.5 * step


def test_plot_mtf_without_a_focal_length_keeps_its_axis():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from KrakenOS.PSFCalc import plot_mtf

    plt.close("all")
    plot_mtf(np.ones((64, 64)), 8.0, 0.55)
    end = float(plt.gca().get_lines()[0].get_xdata()[-1])
    plt.close("all")
    assert end == pytest.approx((8.0 / (0.55 * 1e-3)) / 10.0)
