"""Atmosphere model regressions.

- Mathar computes in every band: above 4.3 um its reference wavenumber was set as `σref` and read
  as `sref` (UnboundLocalError), and the 16-24 um band's cT[3] read -0.954584-9 (= -9.954584) for
  Mathar's -0.954584e-9.
- slalib_refraction computes: it passed `phi` where it had set `_phi` (NameError on every call).
- Observatory.dn_tph takes its own default wavelength uncertainty (dl=None multiplied a float).
"""

import numpy as np
import pytest


def test_mathar_computes_in_every_band():
    from KrakenOS.AstroAtmosphere.refractivityModels import Mathar

    refractivity = [float(np.ravel(Mathar(l))[0]) - 1.0 for l in (2.0, 3.5, 5.0, 10.0, 16.5, 20.0, 23.5)]
    assert all(2.6e-4 < value < 2.9e-4 for value in refractivity), refractivity
    assert np.ptp(refractivity) < 1e-6, refractivity          # dry air barely disperses in the infrared


def test_slalib_refraction_computes():
    from KrakenOS.AstroAtmosphere.misc import slalib_refraction

    refraction = float(np.ravel(slalib_refraction([0.55], 45.0))[0])
    assert refraction == pytest.approx(58.0 / 3600.0, rel=0.05)    # about a minute of arc at 45 degrees


def test_ciddor_uncertainty_takes_its_defaults():
    from KrakenOS.AstroAtmosphere.ciddorModel import Observatory

    uncertainty = float(np.ravel(Observatory().dn_tph(0.55, 288.15, 101325, 0.0, 450))[0])
    assert np.isfinite(uncertainty) and 0.0 < uncertainty < 1e-5
