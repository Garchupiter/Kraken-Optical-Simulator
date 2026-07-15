import math
import time

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


def generate_rays(count):
    side = math.ceil(math.sqrt(count))
    span = np.linspace(-6.0, 6.0, side)
    rays = []
    index = 0
    for x in span:
        for y in span:
            if index >= count:
                break
            rays.append(([float(x), float(y), 0.0], [0.0, 0.0, 1.0], 0.55))
            index += 1
        if index >= count:
            break
    return rays


def extract_minimal_result(system, ray_index):
    xyz = np.asarray(system.XYZ, dtype=float)
    last_lmn = None
    if system.R_LMN:
        last_lmn = np.asarray(system.R_LMN[-1], dtype=float).tolist()

    return {
        "index": ray_index,
        "valid": bool(system.val),
        "last_xyz": xyz[-1].tolist(),
        "last_lmn": last_lmn,
        "top": float(system.TOP),
        "n_hits": int(len(system.XYZ)),
        "surfaces": np.asarray(system.SURFACE).tolist(),
    }


def time_trace_only(rays):
    system = build_simple_system(build=0)
    start = time.perf_counter()
    for origin, direction, wavelength in rays:
        system.Trace(origin, direction, wavelength)
    return time.perf_counter() - start


def time_trace_with_minimal_extract(rays):
    system = build_simple_system(build=0)
    results = []
    start = time.perf_counter()
    for index, (origin, direction, wavelength) in enumerate(rays):
        system.Trace(origin, direction, wavelength)
        results.append(extract_minimal_result(system, index))
    return results, time.perf_counter() - start


def time_trace_with_raykeeper_push(rays):
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    raykeeper = Kos.raykeeper(system)
    start = time.perf_counter()
    for origin, direction, wavelength in rays:
        system.Trace(origin, direction, wavelength)
        raykeeper.push()
    return raykeeper, time.perf_counter() - start


def test_trace_result_collection_costs_are_measurable():
    ray_count = 1000
    rays = generate_rays(ray_count)

    trace_only_time = time_trace_only(rays)
    minimal_results, minimal_time = time_trace_with_minimal_extract(rays)
    raykeeper, raykeeper_time = time_trace_with_raykeeper_push(rays)

    assert len(minimal_results) == ray_count
    assert raykeeper.nrays == ray_count
    assert len(raykeeper.XYZ) == ray_count

    for sample_index in [0, ray_count // 2, ray_count - 1]:
        minimal_last_xyz = minimal_results[sample_index]["last_xyz"]
        raykeeper_last_xyz = np.asarray(raykeeper.XYZ[sample_index], dtype=float)[-1]
        assert np.allclose(minimal_last_xyz, raykeeper_last_xyz, rtol=1e-10, atol=1e-10)

    print(
        "\ntrace collection cost prototype:"
        f"\n  rays={ray_count}"
        f"\n  trace_only={trace_only_time:.6f}s"
        f"\n  trace_plus_minimal_extract={minimal_time:.6f}s"
        f"\n  trace_plus_raykeeper_push={raykeeper_time:.6f}s"
        f"\n  minimal_extract_overhead={minimal_time - trace_only_time:.6f}s"
        f"\n  raykeeper_push_overhead={raykeeper_time - trace_only_time:.6f}s"
    )


def test_raykeeper_ingestion_scales_near_linearly():
    import KrakenOS as Kos
    import tracemalloc

    def ingest(count):
        system = build_simple_system(build=0)
        keeper = Kos.raykeeper(system)
        payload = generate_rays(count)
        tracemalloc.start()
        start = time.perf_counter()
        for origin, direction, wavelength in payload:
            system.Trace(origin, direction, wavelength)
            keeper.push()
        elapsed = time.perf_counter() - start
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        retained = sum(
            np.asarray(item).nbytes if isinstance(item, np.ndarray) else 0
            for item in keeper.XYZ
        )
        return elapsed, peak, current, retained, keeper.nrays

    t1, p1, _, r1, n1 = ingest(500)
    t2, p2, _, r2, n2 = ingest(1000)

    assert n1 == 500
    assert n2 == 1000
    assert t2 / max(t1, 1e-9) < 2.8
    assert p2 / max(p1, 1) < 2.8
    assert r2 / max(r1, 1) < 2.8

    print(
        "\nraykeeper ingestion scaling:"
        f"\n  n=500 time={t1:.6f}s peak={p1} retained_xyz={r1}"
        f"\n  n=1000 time={t2:.6f}s peak={p2} retained_xyz={r2}"
    )
