"""Benchmark warm parallel sequential tracing with worker-local KrakenOS systems.

This is an exploratory benchmark, not a correctness test. It measures:

- sequential tracing with one already-built ``system(build=0)``;
- parallel total time, including process startup and worker initialization;
- parallel warm trace time, after workers have initialized their own systems.

Run from the repository root:

    python tools/benchmark_parallel_trace.py

Optional:

    python tools/benchmark_parallel_trace.py --rays 100 1000 5000 --workers 1 2 4 8
"""

import argparse
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


def trace_sequential(rays):
    system = build_simple_system(build=0)
    start = time.perf_counter()
    results = trace_batch_with_system(system, rays)
    elapsed = time.perf_counter() - start
    return sorted(results, key=lambda item: item["index"]), elapsed


def trace_parallel(rays, workers, batch_size):
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
    if len(sequential) != len(parallel):
        raise AssertionError("result length mismatch")
    for seq, par in zip(sequential, parallel):
        if seq["index"] != par["index"] or seq["val"] != par["val"]:
            raise AssertionError("result identity mismatch")
        if seq["SURFACE"] != par["SURFACE"]:
            raise AssertionError("surface sequence mismatch")
        if not np.allclose(
            np.asarray(seq["XYZ"], dtype=float)[-1],
            np.asarray(par["XYZ"], dtype=float)[-1],
            rtol=1e-10,
            atol=1e-10,
        ):
            raise AssertionError("last XYZ mismatch")
        if not np.isclose(
            np.asarray(seq["TOP"], dtype=float),
            np.asarray(par["TOP"], dtype=float),
            rtol=1e-10,
            atol=1e-10,
        ):
            raise AssertionError("TOP mismatch")

        seq_lmn = result_last_lmn(seq)
        par_lmn = result_last_lmn(par)
        if seq_lmn is None or par_lmn is None:
            if seq_lmn is not None or par_lmn is not None:
                raise AssertionError("last LMN mismatch")
        elif not np.allclose(seq_lmn, par_lmn, rtol=1e-10, atol=1e-10):
            raise AssertionError("last LMN mismatch")


def default_worker_counts():
    cpu_count = os.cpu_count() or 2
    candidates = [1, 2, 4, 8, 16, cpu_count]
    return sorted({workers for workers in candidates if 1 <= workers <= cpu_count})


def run_benchmark(ray_counts, worker_counts):
    print("KrakenOS parallel trace benchmark")
    print(f"cpu_count={os.cpu_count()}")
    print(
        "rays  workers  batch_size  sequential_warm  parallel_total  "
        "parallel_warm_trace  total_speedup  warm_speedup"
    )

    for ray_count in ray_counts:
        rays = generate_rays(ray_count)
        sequential, sequential_time = trace_sequential(rays)

        for workers in worker_counts:
            if workers > (os.cpu_count() or workers):
                continue
            batch_size = max(1, math.ceil(ray_count / workers))
            parallel, parallel_total_time, parallel_trace_time = trace_parallel(
                rays,
                workers=workers,
                batch_size=batch_size,
            )
            assert_results_match(sequential, parallel)

            total_speedup = sequential_time / parallel_total_time if parallel_total_time else float("inf")
            warm_speedup = sequential_time / parallel_trace_time if parallel_trace_time else float("inf")
            print(
                f"{ray_count:4d}  {workers:7d}  {batch_size:10d}  "
                f"{sequential_time:15.6f}s  {parallel_total_time:14.6f}s  "
                f"{parallel_trace_time:19.6f}s  {total_speedup:13.3f}x  "
                f"{warm_speedup:12.3f}x"
            )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rays", nargs="+", type=int, default=[100, 1000])
    parser.add_argument("--workers", nargs="+", type=int, default=default_worker_counts())
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_benchmark(args.rays, args.workers)
