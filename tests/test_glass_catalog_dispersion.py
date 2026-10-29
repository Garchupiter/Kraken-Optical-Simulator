"""Glass catalog and dispersion regressions.

- Formula 3 (Herzberger) is L = 1 / (lambda^2 - 0.028); it was computed with 2 * lambda.
- Every glass of the shipped catalogs gives an index: glasses with no internal-transmittance
  table (the Heraeus fused silicas, the liquids) and the three Zeon plastics, whose NM row is
  shorter, raised instead.
- A section's data wrapped onto following lines (Nikon's CD and TD) is read whole; it was cut to
  its first line, and 24 of Nikon's 137 glasses raised IndexError reaching for a later coefficient.
- The last glass of the catalogs is loaded with its data.
- Setup.Load leaves the caller's list as it was and reads UTILIDADES.AGF once.
- An unknown glass raises a ValueError (it raised a str, a TypeError in Python 3).
"""

import numpy as np
import pytest

D_LINE_UM = 0.5875618


@pytest.fixture(scope="module")
def setup():
    import KrakenOS as Kos

    return Kos.Setup()


def _index(setup, name, wave):
    from KrakenOS.Physics import n_wave_dispersion

    n, _alpha = n_wave_dispersion(setup, name, wave)
    return float(n)


def _first_rows(setup):
    """Each name's first row -- the one a lookup by name finds (some names are in two catalogs)."""
    rows = {}
    for r, name in enumerate(setup.NAMES):
        rows.setdefault(str(name), r)
    return rows


def test_herzberger_formula_takes_the_wavelength_squared(setup):
    """Herzberger glasses that state N(d) inside their range give it back."""
    checked = []
    for name, r in _first_rows(setup).items():
        nm, ld = setup.NM[r], setup.LD[r]
        if len(nm) < 3 or int(nm[0]) != 3 or len(ld) < 2:
            continue
        if not (ld[0] <= D_LINE_UM <= ld[1]) or nm[2] <= 1.0:
            continue
        checked.append((name, abs(_index(setup, name, D_LINE_UM) - nm[2])))
    assert len(checked) >= 5, checked
    worst = max(checked, key=lambda item: item[1])
    assert worst[1] < 1e-4, worst


def test_herzberger_index_of_amtir1_at_10_um(setup):
    assert _index(setup, "AMTIR1", 10.0) == pytest.approx(2.4975, abs=5e-4)


def test_wrapped_catalog_data_is_read_whole(setup):
    """Nikon's glasses (formula 13) wrap their coefficients onto following lines."""
    checked = []
    for name, r in _first_rows(setup).items():
        nm = setup.NM[r]
        if not name.startswith("J-") or len(nm) < 3 or int(nm[0]) != 13:
            continue
        checked.append((name, abs(_index(setup, name, D_LINE_UM) - nm[2])))
    assert len(checked) >= 100, len(checked)
    worst = max(checked, key=lambda item: item[1])
    assert worst[1] < 1e-5, worst


def test_every_catalog_glass_gives_an_index(setup):
    failed = []
    for r, name in enumerate(setup.NAMES):
        ld = setup.LD[r]
        wave = float(np.mean(ld[:2])) if len(ld) >= 2 and ld[1] > ld[0] > 0 else D_LINE_UM
        try:
            n = _index(setup, name, wave)
        except Exception as exc:
            failed.append(f"{name}: {type(exc).__name__}")
            continue
        if not np.isfinite(n):
            failed.append(f"{name}: {n}")
    assert not failed, f"{len(failed)} of {len(setup.NAMES)}: {failed[:12]}"


def test_the_last_glass_of_the_catalogs_is_loaded(setup):
    assert len(setup.NM) == len(setup.NAMES) == len(setup.CD) == len(setup.IT)
    assert np.isfinite(_index(setup, setup.NAMES[-1], D_LINE_UM))


def test_load_leaves_the_callers_list_alone(setup):
    from KrakenOS.SetupClass import rute

    catalogs = [rute + "SCHOTT.AGF"]
    setup.Load(catalogs)
    assert catalogs == [rute + "SCHOTT.AGF"]
    assert sum(1 for name in setup.CAT if name.endswith("UTILIDADES.AGF")) == 1
    setup.Load(list(setup.GlassCat))
    assert sum(1 for name in setup.CAT if name.endswith("UTILIDADES.AGF")) == 1


def test_an_unknown_glass_raises_a_value_error(setup):
    with pytest.raises(ValueError):
        _index(setup, "NO_SUCH_GLASS", D_LINE_UM)
