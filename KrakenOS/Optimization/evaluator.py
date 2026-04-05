from __future__ import annotations

import copy
import io
from contextlib import redirect_stderr, redirect_stdout
from typing import Iterable

import numpy as np

import KrakenOS as Kos

from .merit import MeritFunction, MeritResult, OperandResult
from .operands import (
    EffectiveFocalLengthOperand,
    InvalidTracePenaltyOperand,
    SpotRMSOperand,
    ThicknessPenaltyOperand,
    WavefrontRMSOperand,
)
from .variables import OpticalVariable


class MeritEvaluator:
    """Evaluate merit functions on Kraken surface prescriptions."""

    def __init__(
        self,
        surfaces: Iterable,
        setup=None,
        merit_function: MeritFunction | None = None,
    ) -> None:
        self.base_surfaces = [copy.deepcopy(surface) for surface in surfaces]
        self.setup = setup if setup is not None else Kos.Setup()
        self.merit_function = merit_function if merit_function is not None else MeritFunction([])

    def apply_variables(self, variables: list[OpticalVariable], x: Iterable[float]):
        surfaces = [copy.deepcopy(surface) for surface in self.base_surfaces]
        for variable, value in zip(variables, x):
            if not variable.enabled:
                continue
            if variable.surface_index < 0 or variable.surface_index >= len(surfaces):
                raise IndexError(f"Surface index out of range: {variable.surface_index}")
            setattr(surfaces[variable.surface_index], variable.parameter, float(value))
        return surfaces

    def evaluate(self, variables: list[OpticalVariable], x: Iterable[float]) -> MeritResult:
        surfaces = self.apply_variables(variables, x)
        prescription = [self._surface_to_dict(surface) for surface in surfaces]

        try:
            sink = io.StringIO()
            with redirect_stdout(sink), redirect_stderr(sink):
                system = Kos.system(surfaces, self.setup, build=1)
        except Exception as exc:
            return MeritResult(
                total=self.merit_function.invalid_penalty,
                operands=[],
                valid=False,
                message=f"System build failed: {exc}",
                prescription=prescription,
            )

        results: list[OperandResult] = []
        total = 0.0
        try:
            for operand in self.merit_function.operands:
                if not getattr(operand, "enabled", True):
                    continue
                result = self._evaluate_operand(system, surfaces, operand)
                results.append(result)
                total += result.weighted
        except Exception as exc:
            return MeritResult(
                total=self.merit_function.invalid_penalty,
                operands=results,
                valid=False,
                message=f"Operand evaluation failed: {exc}",
                prescription=prescription,
            )

        return MeritResult(
            total=total,
            operands=results,
            valid=True,
            message="ok",
            prescription=prescription,
        )

    def _evaluate_operand(self, system, surfaces, operand) -> OperandResult:
        if isinstance(operand, SpotRMSOperand):
            value, metadata = self._spot_rms(system, operand)
            residual = value - operand.target
            return OperandResult(
                name=operand.resolved_name(),
                value=value,
                target=operand.target,
                residual=residual,
                weighted=operand.weight * residual * residual,
                metadata=metadata,
            )

        if isinstance(operand, EffectiveFocalLengthOperand):
            _, _, _, _, _, _, _, effl, _, _, _, _, _ = system.Parax(operand.wavelength)
            value = float(effl)
            residual = value - operand.target
            return OperandResult(
                name=operand.resolved_name(),
                value=value,
                target=operand.target,
                residual=residual,
                weighted=operand.weight * residual * residual,
                metadata={},
            )

        if isinstance(operand, WavefrontRMSOperand):
            value, metadata = self._wavefront_rms(system, operand)
            residual = value - operand.target
            return OperandResult(
                name=operand.resolved_name(),
                value=value,
                target=operand.target,
                residual=residual,
                weighted=operand.weight * residual * residual,
                metadata=metadata,
            )

        if isinstance(operand, ThicknessPenaltyOperand):
            thickness = float(surfaces[operand.surface_index].Thickness)
            penalty = 0.0
            if thickness < operand.minimum:
                penalty += (operand.minimum - thickness) ** 2
            if operand.maximum is not None and thickness > operand.maximum:
                penalty += (thickness - operand.maximum) ** 2
            return OperandResult(
                name=operand.resolved_name(),
                value=thickness,
                target=operand.target,
                residual=penalty,
                weighted=operand.weight * penalty,
                metadata={"minimum": operand.minimum, "maximum": operand.maximum},
            )

        if isinstance(operand, InvalidTracePenaltyOperand):
            return OperandResult(
                name=operand.resolved_name(),
                value=0.0,
                target=operand.target,
                residual=0.0,
                weighted=0.0,
                metadata={"penalty": operand.penalty},
            )

        raise TypeError(f"Unsupported operand type: {type(operand)!r}")

    def _spot_rms(self, system, operand: SpotRMSOperand):
        rays = Kos.raykeeper(system)
        max_radius = max((float(surface.Diameter) / 2.0 for surface in system.SDT), default=1.0)
        span = max_radius * max(min(operand.ray_height_factor, 1.5), 0.05)
        ray_count = max(1, int(operand.ray_count))
        heights = [0.0] if ray_count == 1 else np.linspace(-span, span, ray_count)

        for y0 in heights:
            system.Trace([0.0, float(y0), 0.0], [0.0, 0.0, 1.0], operand.wavelength)
            rays.push()

        X, Y, Z, L, M, N = rays.pick(operand.surface_index)
        if X.size == 0:
            raise RuntimeError("No valid rays reached the requested surface")
        rms, cen_x, cen_y = Kos.RMS(X, Y, Z, L, M, N)
        return float(rms), {
            "centroid_x": float(cen_x),
            "centroid_y": float(cen_y),
            "surface_index": operand.surface_index,
            "ray_count": ray_count,
        }

    def _wavefront_rms(self, system, operand: WavefrontRMSOperand):
        pupil = Kos.PupilCalc(
            system,
            operand.surface_index,
            operand.wavelength,
            operand.aperture_type,
            operand.aperture_value,
        )
        pupil.Samp = int(operand.sample_size)
        pupil.Ptype = operand.pattern
        pupil.FieldType = operand.field_type
        pupil.FieldX = float(operand.field_x)
        pupil.FieldY = float(operand.field_y)

        x, y, phase, p2v = Kos.Phase(pupil)
        phase = np.asarray(phase, dtype=float)
        finite = np.isfinite(phase)
        if not np.any(finite):
            raise RuntimeError("Wavefront phase evaluation returned no finite samples")

        finite_phase = phase[finite]
        centered = finite_phase - np.mean(finite_phase)
        rms = float(np.sqrt(np.mean(centered * centered)))
        return rms, {
            "surface_index": operand.surface_index,
            "sample_size": int(operand.sample_size),
            "peak_to_valley": float(p2v),
            "sample_count": int(finite_phase.size),
        }

    @staticmethod
    def _surface_to_dict(surface) -> dict:
        return {
            "Name": getattr(surface, "Name", ""),
            "Rc": float(getattr(surface, "Rc", 0.0)),
            "Thickness": float(getattr(surface, "Thickness", 0.0)),
            "Diameter": float(getattr(surface, "Diameter", 0.0)),
            "Glass": getattr(surface, "Glass", "AIR"),
            "Thin_Lens": float(getattr(surface, "Thin_Lens", 0.0)),
        }
