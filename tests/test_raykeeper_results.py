import tracemalloc

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


def trace_valid_ray(system):
    system.Trace([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.55)


def trace_invalid_ray(system):
    system.Trace([1000.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.55)


def assert_array_lists_equal(left, right):
    assert len(left) == len(right)
    for left_item, right_item in zip(left, right):
        np.testing.assert_equal(left_item, right_item)


def assert_raykeepers_match(left, right):
    assert left.nrays == right.nrays
    np.testing.assert_equal(left.vld, right.vld)

    list_attrs = [
        "RayWave",
        "CC",
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
        "valid_SURFACE",
        "valid_NAME",
        "valid_GLASS",
        "valid_XYZ",
        "valid_OST_XYZ",
        "valid_LMN",
        "valid_R_LMN",
        "valid_TOP",
        "valid_ALPHA",
        "invalid_SURFACE",
        "invalid_NAME",
        "invalid_GLASS",
        "invalid_XYZ",
        "invalid_OST_XYZ",
        "invalid_LMN",
        "invalid_R_LMN",
        "invalid_TOP",
        "invalid_ALPHA",
    ]
    for attr in list_attrs:
        assert_array_lists_equal(getattr(left, attr), getattr(right, attr))


def test_push_result_matches_push_for_valid_trace():
    import KrakenOS as Kos

    push_system = build_simple_system(build=0)
    result_system = build_simple_system(build=0)
    push_rays = Kos.raykeeper(push_system)
    result_rays = Kos.raykeeper(result_system)

    trace_valid_ray(push_system)
    push_rays.push()

    trace_valid_ray(result_system)
    result = Kos.extract_ray_result(result_system)
    result_rays.push_result(result)

    assert_raykeepers_match(push_rays, result_rays)


def test_push_result_matches_push_for_invalid_trace():
    import KrakenOS as Kos

    push_system = build_simple_system(build=0)
    result_system = build_simple_system(build=0)
    push_rays = Kos.raykeeper(push_system)
    result_rays = Kos.raykeeper(result_system)

    trace_invalid_ray(push_system)
    push_rays.push()

    trace_invalid_ray(result_system)
    result = Kos.extract_ray_result(result_system)
    result_rays.push_result(result)

    assert_raykeepers_match(push_rays, result_rays)


def test_extend_results_stores_multiple_extracted_results():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    results = []

    trace_valid_ray(system)
    results.append(Kos.extract_ray_result(system, copy=True))

    trace_invalid_ray(system)
    results.append(Kos.extract_ray_result(system, copy=True))

    rays.extend_results(results)

    assert rays.nrays == 2
    assert rays.vld.tolist() == [1.0]
    assert len(rays.valid_XYZ) == 1
    assert len(rays.invalid_XYZ) == 1
    assert len(rays.XYZ) == 2


def test_mixed_valid_invalid_order_and_buckets():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)

    trace_valid_ray(system)
    rays.push()
    trace_invalid_ray(system)
    rays.push()
    trace_valid_ray(system)
    rays.push()

    assert rays.nrays == 3
    assert rays.vld.tolist() == [1.0, 1.0]
    assert len(rays.valid_XYZ) == 2
    assert len(rays.invalid_XYZ) == 1
    assert len(rays.XYZ) == 3
    assert len(list(rays.valid())) == 2
    # universal order preserves push order; buckets filter by validity
    assert np.allclose(rays.XYZ[0].astype(float)[-1], rays.valid_XYZ[0].astype(float)[-1])
    assert np.allclose(rays.XYZ[2].astype(float)[-1], rays.valid_XYZ[1].astype(float)[-1])
    assert np.allclose(rays.XYZ[1].astype(float), rays.invalid_XYZ[0].astype(float))


def test_ost_xyz_and_alpha_consistent_across_buckets():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    trace_valid_ray(system)
    rays.push()

    assert len(rays.OST_XYZ) == 1
    assert len(rays.valid_OST_XYZ) == 1
    np.testing.assert_equal(rays.OST_XYZ[0], rays.valid_OST_XYZ[0])
    np.testing.assert_equal(rays.ALPHA[0], rays.valid_ALPHA[0])


def test_extract_snapshot_survives_next_trace():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)

    trace_valid_ray(system)
    snapshot = Kos.extract_ray_result(system, copy=True)
    first_xyz = np.asarray(snapshot["XYZ"], dtype=float).copy()

    trace_invalid_ray(system)
    rays.push_result(snapshot)

    assert rays.nrays == 1
    assert rays.vld.tolist() == [1.0]
    assert np.allclose(np.asarray(rays.XYZ[0], dtype=float), first_xyz)


def test_clean_resets_state():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    trace_valid_ray(system)
    rays.push()
    rays.pick(-1)

    rays.clean()

    assert rays.nrays == 0
    assert rays.nelements == 0
    assert rays.vld.size == 0
    assert len(rays.XYZ) == 0
    assert len(rays.valid_XYZ) == 0
    assert len(rays.invalid_XYZ) == 0
    assert getattr(rays, "valid_vld").size == 0
    assert getattr(rays, "invalid_vld").size == 0


def test_shared_storage_reuses_payload_refs():
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    trace_valid_ray(system)
    rays.push()

    assert rays.XYZ[0] is rays.valid_XYZ[0]
    assert rays.SURFACE[0] is rays.valid_SURFACE[0]
    assert rays.OST_XYZ[0] is rays.valid_OST_XYZ[0]


def test_raykeeper_push_memory_grows_near_linearly():
    import KrakenOS as Kos
    import time

    def measure(count):
        system = build_simple_system(build=0)
        rays = Kos.raykeeper(system)
        tracemalloc.start()
        start = time.perf_counter()
        for _ in range(count):
            trace_valid_ray(system)
            rays.push()
        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return elapsed, peak, rays.nrays

    time_small, peak_small, n_small = measure(400)
    time_large, peak_large, n_large = measure(800)

    assert n_small == 400
    assert n_large == 800
    # Contiguous np.append would grow roughly quadratically; require near-linear.
    assert time_large / max(time_small, 1e-9) < 2.8
    assert peak_large / max(peak_small, 1) < 2.8
