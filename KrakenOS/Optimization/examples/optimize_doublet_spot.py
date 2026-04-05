#!/usr/bin/env python3

from __future__ import annotations

import json
import os
import sys
import ctypes
from pathlib import Path

from scipy.optimize import differential_evolution

import KrakenOS as Kos

from KrakenOS.Optimization import (
    EffectiveFocalLengthOperand,
    MeritEvaluator,
    MeritFunction,
    OpticalVariable,
    SpotRMSOperand,
    ThicknessPenaltyOperand,
    WavefrontRMSOperand,
)
from KrakenOS.Optimization.adapters.pygmo2_adapter import Pygmo2MeritProblem

SCIPY_MAXITER = 4
SCIPY_POPSIZE = 4
PYGMO_POPSIZE = 8
PYGMO_GENERATIONS = 8
PYGMO_VERBOSITY = 2


def build_doublet():
    p_obj = Kos.surf()
    p_obj.Name = "Object"
    p_obj.Rc = 0.0
    p_obj.Thickness = 10.0
    p_obj.Glass = "AIR"
    p_obj.Diameter = 30.0

    l1a = Kos.surf()
    l1a.Name = "Crown Front"
    l1a.Rc = 92.8470657
    l1a.Thickness = 6.0
    l1a.Glass = "BK7"
    l1a.Diameter = 30.0

    l1b = Kos.surf()
    l1b.Name = "Flint Front"
    l1b.Rc = -30.7160867
    l1b.Thickness = 3.0
    l1b.Glass = "F2"
    l1b.Diameter = 30.0

    l1c = Kos.surf()
    l1c.Name = "Flint Back"
    l1c.Rc = -78.19730726
    l1c.Thickness = 97.37604743
    l1c.Glass = "AIR"
    l1c.Diameter = 30.0

    p_img = Kos.surf()
    p_img.Name = "Image"
    p_img.Rc = 0.0
    p_img.Thickness = 0.0
    p_img.Glass = "AIR"
    p_img.Diameter = 3.0

    return [p_obj, l1a, l1b, l1c, p_img]


def build_problem():
    surfaces = build_doublet()
    setup = Kos.Setup()

    variables = [
        OpticalVariable(1, "Rc", 60.0, 120.0, name="Crown Front Rc"),
        OpticalVariable(2, "Rc", -60.0, -10.0, name="Flint Front Rc"),
        OpticalVariable(3, "Rc", -120.0, -40.0, name="Flint Back Rc"),
        OpticalVariable(3, "Thickness", 60.0, 120.0, name="Back Focus"),
    ]

    merit = MeritFunction(
        operands=[
            SpotRMSOperand(
                name="Spot RMS",
                weight=1.0,
                target=0.0,
                surface_index=-1,
                wavelength=0.55,
                ray_count=11,
            ),
            EffectiveFocalLengthOperand(
                name="EFFL target",
                weight=1e-4,
                target=100.0,
                wavelength=0.55,
            ),
            WavefrontRMSOperand(
                name="Wavefront RMS",
                weight=1e-2,
                target=0.0,
                surface_index=1,
                wavelength=0.55,
                aperture_type="EPD",
                aperture_value=30.0,
                sample_size=9,
            ),
            ThicknessPenaltyOperand(
                name="Back focus bounds",
                weight=1.0,
                surface_index=3,
                minimum=10.0,
                maximum=150.0,
            ),
        ]
    )

    evaluator = MeritEvaluator(surfaces, setup=setup, merit_function=merit)
    # Start from a deliberately perturbed prescription so the optimization
    # run demonstrates a meaningful improvement in the recorded champion.
    x0 = [110.0, -20.0, -60.0, 80.0]
    return evaluator, variables, x0


def result_to_dict(result, x):
    return {
        "x": [float(v) for v in x],
        "total": float(result.total),
        "valid": bool(result.valid),
        "message": result.message,
        "operands": [
            {
                "name": operand.name,
                "value": float(operand.value),
                "target": float(operand.target),
                "residual": float(operand.residual),
                "weighted": float(operand.weighted),
                "metadata": operand.metadata,
            }
            for operand in result.operands
        ],
        "prescription": result.prescription,
    }


def run_scipy(evaluator, variables, x0):
    def objective(x):
        return evaluator.evaluate(variables, x).total

    bounds = [(variable.lower_bound, variable.upper_bound) for variable in variables]
    de_result = differential_evolution(
        objective,
        bounds=bounds,
        maxiter=SCIPY_MAXITER,
        popsize=SCIPY_POPSIZE,
        polish=False,
        seed=42,
    )
    merit_result = evaluator.evaluate(variables, de_result.x)
    return {
        "backend": "scipy.differential_evolution",
        "status": de_result.message,
        "nit": int(de_result.nit),
        "nfev": int(de_result.nfev),
        "fun": float(de_result.fun),
        "initial_total": float(evaluator.evaluate(variables, x0).total),
        "result": result_to_dict(merit_result, de_result.x),
    }


def import_local_pygmo():
    pagmo_lib = Path(os.path.expanduser("~/Projects/pagmo2/_install/lib64/libpagmo.so"))
    if pagmo_lib.exists():
        try:
            ctypes.CDLL(str(pagmo_lib), mode=ctypes.RTLD_GLOBAL)
        except OSError:
            pass

    try:
        import pygmo as pg  # type: ignore

        return pg, None
    except Exception as exc:
        direct_error = exc

    repo_root = Path(os.path.expanduser("~/Projects/pygmo2"))
    if repo_root.exists():
        sys.path.insert(0, str(repo_root))
        try:
            import pygmo as pg  # type: ignore

            return pg, None
        except Exception as exc:
            return None, f"{direct_error}; local repo import failed: {exc}"

    return None, str(direct_error)


def run_pygmo(evaluator, variables, x0):
    pg, error = import_local_pygmo()
    if pg is None:
        return {
            "backend": "pygmo",
            "available": False,
            "error": error,
        }

    udp = Pygmo2MeritProblem(evaluator=evaluator, variables=variables)
    problem = pg.problem(udp)
    population = pg.population(problem, size=PYGMO_POPSIZE, seed=42)
    population.push_back(x0)
    algorithm = pg.algorithm(pg.de(gen=PYGMO_GENERATIONS, seed=42))
    algorithm.set_verbosity(PYGMO_VERBOSITY)
    evolved = algorithm.evolve(population)
    log = [
        {
            "gen": int(gen),
            "fevals": int(fevals),
            "best": float(best),
            "dx": float(dx),
            "df": float(df),
        }
        for gen, fevals, best, dx, df in algorithm.extract(pg.de).get_log()
    ]
    x_best = evolved.champion_x
    merit_result = evaluator.evaluate(variables, x_best)
    return {
        "backend": "pygmo.de",
        "available": True,
        "population_size": int(len(evolved.get_x())),
        "generations": PYGMO_GENERATIONS,
        "initial_total": float(evaluator.evaluate(variables, x0).total),
        "log": log,
        "result": result_to_dict(merit_result, x_best),
    }


def main():
    evaluator, variables, x0 = build_problem()
    initial = evaluator.evaluate(variables, x0)

    print("Initial merit:", initial.total)
    for operand in initial.operands:
        print(f"{operand.name}: value={operand.value:.6g} weighted={operand.weighted:.6g}")

    problem = Pygmo2MeritProblem(evaluator=evaluator, variables=variables)
    print("Bounds:", problem.get_bounds())

    outputs = {
        "initial": result_to_dict(initial, x0),
        "scipy": run_scipy(evaluator, variables, x0),
        "pygmo": run_pygmo(evaluator, variables, x0),
    }

    if outputs["pygmo"].get("available"):
        print(
            "pygmo improvement:",
            outputs["pygmo"]["initial_total"],
            "->",
            outputs["pygmo"]["result"]["total"],
        )
    print(
        "scipy improvement:",
        outputs["scipy"]["initial_total"],
        "->",
        outputs["scipy"]["result"]["total"],
    )

    out_path = Path("doublet_merit_result.json")
    out_path.write_text(json.dumps(outputs, indent=2) + "\n", encoding="utf-8")
    print("Wrote", out_path)


if __name__ == "__main__":
    main()
