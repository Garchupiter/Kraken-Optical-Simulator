import math
import os
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np


_WORKER_SYSTEM = None


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


def init_worker_system():
    global _WORKER_SYSTEM
    _WORKER_SYSTEM = build_simple_system(build=0)


def generate_rays(count):
    side = math.ceil(math.sqrt(count))
    span = np.linspace(-6.0, 6.0, side)
    rays = []
    index = 0
    for x in span:
        for y in span:
            if index >= count:
                break
            rays.append(
                {
                    "index": index,
                    "origin": [float(x), float(y), 0.0],
                    "direction": [0.0, 0.0, 1.0],
                    "wavelength": 0.55,
                }
            )
            index += 1
        if index >= count:
            break
    return rays


def make_warmup_batches(workers):
    return [
        [
            {
                "index": -worker_index - 1,
                "origin": [0.0, 0.0, 0.0],
                "direction": [0.0, 0.0, 1.0],
                "wavelength": 0.55,
            }
        ]
        for worker_index in range(workers)
    ]


def trace_batch(batch):
    system = build_simple_system(build=0)
    return trace_batch_with_system(system, batch)


def trace_batch_with_system(system, batch):
    import KrakenOS as Kos

    results = []
    for ray in batch:
        system.Trace(ray["origin"], ray["direction"], ray["wavelength"])
        result = Kos.extract_ray_result(system, copy=True)
        result["index"] = ray["index"]
        results.append(result)
    return results


def trace_batch_with_worker_system(batch):
    if _WORKER_SYSTEM is None:
        raise RuntimeError("Worker system was not initialized")
    return trace_batch_with_system(_WORKER_SYSTEM, batch)


def chunked(items, chunk_size):
    return [items[index : index + chunk_size] for index in range(0, len(items), chunk_size)]


def worker_counts_to_test():
    return [min(2, os.cpu_count() or 2)]


def trace_sequential(rays):
    system = build_simple_system(build=0)
    start = time.perf_counter()
    results = trace_batch_with_system(system, rays)
    elapsed = time.perf_counter() - start
    return sorted(results, key=lambda item: item["index"]), elapsed


def trace_parallel(rays, workers=2, batch_size=25):
    batches = chunked(rays, batch_size)
    total_start = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=get_context("spawn"),
        initializer=init_worker_system,
    ) as pool:
        list(pool.map(trace_batch_with_worker_system, make_warmup_batches(workers)))
        trace_start = time.perf_counter()
        grouped_results = list(pool.map(trace_batch_with_worker_system, batches))
        trace_elapsed = time.perf_counter() - trace_start
    total_elapsed = time.perf_counter() - total_start
    results = [result for group in grouped_results for result in group]
    return sorted(results, key=lambda item: item["index"]), total_elapsed, trace_elapsed


def result_last_lmn(result):
    if not result["R_LMN"]:
        return None
    candidate = np.asarray(result["R_LMN"][-1], dtype=float)
    if not candidate.size:
        return None
    return candidate


def assert_results_match(sequential, parallel):
    assert len(sequential) == len(parallel)
    for seq, par in zip(sequential, parallel):
        assert seq["index"] == par["index"]
        assert seq["val"] == par["val"]
        assert len(seq["XYZ"]) == len(par["XYZ"])
        assert seq["SURFACE"] == par["SURFACE"]
        assert np.allclose(np.asarray(seq["XYZ"], dtype=float)[-1], np.asarray(par["XYZ"], dtype=float)[-1], rtol=1e-10, atol=1e-10)
        assert np.isclose(np.asarray(seq["TOP"], dtype=float), np.asarray(par["TOP"], dtype=float), rtol=1e-10, atol=1e-10)

        seq_lmn = result_last_lmn(seq)
        par_lmn = result_last_lmn(par)
        if seq_lmn is None or par_lmn is None:
            assert seq_lmn is None and par_lmn is None
        else:
            assert np.allclose(seq_lmn, par_lmn, rtol=1e-10, atol=1e-10)


def build_raykeeper_from_results(results):
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    rays.extend_results(sorted(results, key=lambda item: item["index"]))
    return rays


def build_raykeeper_from_batches(grouped_results):
    """Ingest batches as they arrive to avoid retaining one giant result list."""
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    for batch in grouped_results:
        rays.extend_results(sorted(batch, key=lambda item: item["index"]))
    return rays


def build_classic_raykeeper(rays_to_trace):
    import KrakenOS as Kos

    system = build_simple_system(build=0)
    rays = Kos.raykeeper(system)
    for ray in rays_to_trace:
        system.Trace(ray["origin"], ray["direction"], ray["wavelength"])
        rays.push()
    return rays


def assert_raykeepers_match(reference, candidate):
    assert reference.nrays == candidate.nrays
    np.testing.assert_equal(reference.vld, candidate.vld)
    assert len(reference.XYZ) == len(candidate.XYZ)
    assert len(reference.R_LMN) == len(candidate.R_LMN)
    assert len(reference.CC) == len(candidate.CC)
    for ref_xyz, cand_xyz in zip(reference.XYZ, candidate.XYZ):
        assert np.allclose(ref_xyz.astype(float), cand_xyz.astype(float), rtol=1e-10, atol=1e-10)
    for ref_lmn, cand_lmn in zip(reference.R_LMN, candidate.R_LMN):
        assert np.allclose(ref_lmn.astype(float), cand_lmn.astype(float), rtol=1e-10, atol=1e-10)
    for ref_hits, cand_hits in zip(reference.CC, candidate.CC):
        assert np.allclose(ref_hits.astype(float), cand_hits.astype(float), rtol=1e-10, atol=1e-10)


def test_parallel_sequential_trace_matches_serial_results():
    ray_count = 100
    rays = generate_rays(ray_count)

    sequential, sequential_time = trace_sequential(rays)
    classic_raykeeper = build_classic_raykeeper(rays)

    print(
        "\nparallel trace prototype:"
        f"\n  rays={ray_count}"
        f"\n  cpu_count={os.cpu_count()}"
        f"\n  sequential_warm={sequential_time:.6f}s"
        "\n  workers  batch_size  parallel_total  parallel_warm_trace  total_speedup  warm_speedup"
    )

    for workers in worker_counts_to_test():
        batch_size = max(1, math.ceil(ray_count / workers))
        batches = chunked(rays, batch_size)
        total_start = time.perf_counter()
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=get_context("spawn"),
            initializer=init_worker_system,
        ) as pool:
            list(pool.map(trace_batch_with_worker_system, make_warmup_batches(workers)))
            trace_start = time.perf_counter()
            grouped_results = list(pool.map(trace_batch_with_worker_system, batches))
            parallel_trace_time = time.perf_counter() - trace_start
        parallel_total_time = time.perf_counter() - total_start

        parallel = [result for group in grouped_results for result in group]
        parallel = sorted(parallel, key=lambda item: item["index"])
        assert_results_match(sequential, parallel)

        reconstructed_raykeeper = build_raykeeper_from_batches(grouped_results)
        assert_raykeepers_match(classic_raykeeper, reconstructed_raykeeper)

        total_speedup = sequential_time / parallel_total_time if parallel_total_time else float("inf")
        warm_speedup = sequential_time / parallel_trace_time if parallel_trace_time else float("inf")
        print(
            f"  {workers:7d}  {batch_size:10d}  "
            f"{parallel_total_time:14.6f}s  {parallel_trace_time:19.6f}s  "
            f"{total_speedup:13.3f}x  {warm_speedup:12.3f}x"
        )


def test_extract_snapshot_copy_is_independent_after_pickle_roundtrip():
    import pickle

    import KrakenOS as Kos

    system = build_simple_system(build=0)
    system.Trace([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.55)
    snapshot = Kos.extract_ray_result(system, copy=True)
    payload = pickle.dumps(snapshot)
    restored = pickle.loads(payload)

    system.Trace([1000.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.55)
    rays = Kos.raykeeper(system)
    rays.push_result(restored)

    assert rays.vld.tolist() == [1.0]
    assert not np.allclose(np.asarray(restored["XYZ"], dtype=float)[-1], [1000.0, 0.0, 0.0])
    assert not np.allclose(np.asarray(rays.XYZ[0], dtype=float)[-1], [1000.0, 0.0, 0.0])
