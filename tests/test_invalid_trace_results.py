import numpy as np


def build_simple_system(build=0):
    import KrakenOS as Kos

    obj = Kos.surf()
    obj.Glass = "AIR"
    obj.Thickness = 10
    obj.Diameter = 30

    lens_front = Kos.surf()
    lens_front.Rc = 80
    lens_front.Glass = "BK7"
    lens_front.Thickness = 5
    lens_front.Diameter = 30

    lens_back = Kos.surf()
    lens_back.Rc = -80
    lens_back.Glass = "AIR"
    lens_back.Thickness = 20
    lens_back.Diameter = 30

    image = Kos.surf()
    image.Glass = "AIR"
    image.Thickness = 0
    image.Diameter = 30

    return Kos.system([obj, lens_front, lens_back, image], Kos.Setup(), build=build)


def trace_invalid_ray(system):
    system.Trace([1000.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.55)


def extract_safe_minimal_result(system):
    ray_hits = np.asarray(system.ray_SurfHits, dtype=float)
    xyz = np.asarray(system.XYZ, dtype=float)
    last_lmn = None
    if system.R_LMN:
        try:
            candidate = np.asarray(system.R_LMN[-1], dtype=float)
        except ValueError:
            candidate = np.asarray([])
        if candidate.size:
            last_lmn = candidate.tolist()

    try:
        top = float(system.TOP)
    except (TypeError, ValueError):
        top = None

    return {
        "valid": bool(system.val),
        "ray_hits": ray_hits.tolist(),
        "last_xyz": xyz[-1].tolist(),
        "last_lmn": last_lmn,
        "top": top,
        "n_hits": int(len(system.XYZ)),
        "surfaces": np.asarray(system.SURFACE).tolist(),
    }


def test_invalid_trace_current_empty_collect_shape_is_documented():
    system = build_simple_system(build=0)

    trace_invalid_ray(system)

    assert system.val == 0
    assert system.SURFACE == [1]
    assert system.GLASS == ["AIR"]
    assert np.allclose(system.ray_SurfHits, [[1000.0, 0.0, 0.0]])
    assert np.allclose(np.asarray(system.XYZ, dtype=float), [[1000.0, 0.0, 0.0], [1000.0, 0.0, 0.0]])

    # These fields document the current __EmptyCollect/__CollectData behavior.
    # They are odd, but recording them lets a future cleanup prove compatibility
    # or intentionally update the expected invalid-ray contract.
    assert np.asarray(system.WAV).size == 0
    assert np.asarray(system.TOP).size == 0
    assert np.asarray(system.N0[0]).size == 0
    assert np.asarray(system.N1[0]).item() == 0.55
    assert system.TT == 0.0


def test_invalid_trace_can_be_reduced_to_safe_serializable_result():
    system = build_simple_system(build=0)

    trace_invalid_ray(system)
    result = extract_safe_minimal_result(system)

    assert result == {
        "valid": False,
        "ray_hits": [[1000.0, 0.0, 0.0]],
        "last_xyz": [1000.0, 0.0, 0.0],
        "last_lmn": None,
        "top": None,
        "n_hits": 2,
        "surfaces": [1],
    }


def test_raykeeper_push_records_invalid_trace_separately():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)

    trace_invalid_ray(system)
    rays.push()

    assert rays.nrays == 1
    assert rays.vld.size == 0
    assert len(rays.invalid_XYZ) == 1
    assert len(rays.valid_XYZ) == 0
    assert np.allclose(rays.CC[0], [[1000.0, 0.0, 0.0]])
    assert np.allclose(rays.XYZ[0].astype(float), [[1000.0, 0.0, 0.0], [1000.0, 0.0, 0.0]])


def test_invalid_trace_does_not_extend_vld_after_valid_rays():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)

    system.Trace([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.55)
    rays.push()
    trace_invalid_ray(system)
    rays.push()

    assert rays.nrays == 2
    assert rays.vld.tolist() == [1.0]
    assert len(rays.valid_XYZ) == 1
    assert len(rays.invalid_XYZ) == 1
    assert rays.XYZ[0] is rays.valid_XYZ[0]
    assert rays.XYZ[1] is rays.invalid_XYZ[0]
