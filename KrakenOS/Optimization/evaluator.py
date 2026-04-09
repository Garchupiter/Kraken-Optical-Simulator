from __future__ import annotations

import copy
import io
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stderr, redirect_stdout
from typing import Iterable

import numpy as np

import KrakenOS as Kos

from .merit import MeritFunction, MeritResult, OperandResult
from .operands import (
    EffectiveFocalLengthOperand,
    EntrancePupilPositionOperand,
    ExitPupilPositionOperand,
    InvalidTracePenaltyOperand,
    MagnificationOperand,
    MTFAtFrequencyOperand,
    SpotRMSOperand,
    ThicknessPenaltyOperand,
    WavefrontRMSOperand,
)
from .variables import OpticalVariable


_WORKER_SYSTEM_CACHE_SIGNATURE = None
_WORKER_SYSTEM_CACHE_SYSTEM = None


def _surface_signature(surface_specs: list[dict]) -> tuple:
    signature = []
    for spec in surface_specs:
        signature.append(
            (
                str(spec.get("name", "")),
                float(spec.get("rc", 0.0)),
                float(spec.get("thickness", 0.0)),
                float(spec.get("diameter", 0.0)),
                str(spec.get("glass", "AIR")),
                float(spec.get("thin_lens", 0.0)),
                float(spec.get("diff_ord", 0.0)),
                float(spec.get("grating_d", 0.0)),
                float(spec.get("tilt_x", 0.0)),
                float(spec.get("tilt_y", 0.0)),
                float(spec.get("tilt_z", 0.0)),
                float(spec.get("desp_x", 0.0)),
                float(spec.get("desp_y", 0.0)),
                float(spec.get("desp_z", 0.0)),
                float(spec.get("axis_move", 0.0)),
                float(spec.get("drawing", 1.0)),
            )
        )
    return tuple(signature)


def _build_system_from_specs(surface_specs: list[dict]):
    surfaces = []
    for spec in surface_specs:
        surface = Kos.surf()
        surface.Name = str(spec.get("name", ""))
        surface.Rc = float(spec.get("rc", 0.0))
        surface.Thickness = float(spec.get("thickness", 0.0))
        surface.Diameter = float(spec.get("diameter", 0.0))
        surface.Glass = str(spec.get("glass", "AIR"))
        surface.TiltX = float(spec.get("tilt_x", 0.0))
        surface.TiltY = float(spec.get("tilt_y", 0.0))
        surface.TiltZ = float(spec.get("tilt_z", 0.0))
        surface.DespX = float(spec.get("desp_x", 0.0))
        surface.DespY = float(spec.get("desp_y", 0.0))
        surface.DespZ = float(spec.get("desp_z", 0.0))
        surface.AxisMove = float(spec.get("axis_move", 0.0))
        surface.Drawing = float(spec.get("drawing", 1.0))
        thin_lens = float(spec.get("thin_lens", 0.0))
        if abs(thin_lens) > 0.0:
            surface.Thin_Lens = thin_lens
        diff_ord = float(spec.get("diff_ord", 0.0))
        if abs(diff_ord) > 0.0:
            surface.Diff_Ord = diff_ord
            surface.Grating_D = float(spec.get("grating_d", 1.0))
        surfaces.append(surface)
    return Kos.system(surfaces, Kos.Setup(), build=1)


def _build_cached_system_from_specs(surface_specs: list[dict]):
    global _WORKER_SYSTEM_CACHE_SIGNATURE, _WORKER_SYSTEM_CACHE_SYSTEM
    signature = _surface_signature(surface_specs)
    if _WORKER_SYSTEM_CACHE_SYSTEM is None or _WORKER_SYSTEM_CACHE_SIGNATURE != signature:
        _WORKER_SYSTEM_CACHE_SYSTEM = _build_system_from_specs(surface_specs)
        _WORKER_SYSTEM_CACHE_SIGNATURE = signature
    return _WORKER_SYSTEM_CACHE_SYSTEM


def _pick_image_plane_data_static(rays):
    try:
        X, Y, Z, L, M, N = rays.pick(-1, coordinates="local")
        if np.asarray(X).size:
            return X, Y, Z, L, M, N
    except Exception:
        pass
    return rays.pick(-1)


def _trace_mtf_chunk(
    surface_specs: list[dict],
    wavelength: float,
    x_bundle,
    y_bundle,
    z_bundle,
    l_bundle,
    m_bundle,
    n_bundle,
):
    system = _build_cached_system_from_specs(surface_specs)
    rays = Kos.raykeeper(system)
    Kos.TraceLoop(
        np.asarray(x_bundle, dtype=float),
        np.asarray(y_bundle, dtype=float),
        np.asarray(z_bundle, dtype=float),
        np.asarray(l_bundle, dtype=float),
        np.asarray(m_bundle, dtype=float),
        np.asarray(n_bundle, dtype=float),
        float(wavelength),
        rays,
        clean=1,
    )
    X, Y, _Z, _L, _M, _N = _pick_image_plane_data_static(rays)
    return np.asarray(X, dtype=float), np.asarray(Y, dtype=float)


class MeritEvaluator:
    """Evaluate merit functions on Kraken surface prescriptions."""

    def __init__(
        self,
        surfaces: Iterable,
        setup=None,
        merit_function: MeritFunction | None = None,
        mtf_worker_count: int = 1,
    ) -> None:
        self.base_surfaces = [copy.deepcopy(surface) for surface in surfaces]
        self.setup = setup if setup is not None else Kos.Setup()
        self.merit_function = merit_function if merit_function is not None else MeritFunction([])
        self.mtf_worker_count = max(1, int(mtf_worker_count))
        self._mtf_executor: ProcessPoolExecutor | None = None
        self._mtf_executor_workers = 0

    def __del__(self):
        try:
            self._shutdown_mtf_executor()
        except Exception:
            pass

    def __getstate__(self):
        state = dict(self.__dict__)
        state["_mtf_executor"] = None
        state["_mtf_executor_workers"] = 0
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._mtf_executor = None
        self._mtf_executor_workers = 0

    def _shutdown_mtf_executor(self) -> None:
        executor = self._mtf_executor
        self._mtf_executor = None
        self._mtf_executor_workers = 0
        if executor is not None:
            executor.shutdown(wait=False, cancel_futures=True)

    def _ensure_mtf_executor(self, worker_count: int) -> ProcessPoolExecutor | None:
        worker_count = max(1, int(worker_count))
        if worker_count <= 1:
            return None
        if self._mtf_executor is not None and self._mtf_executor_workers == worker_count:
            return self._mtf_executor
        self._shutdown_mtf_executor()
        self._mtf_executor = ProcessPoolExecutor(max_workers=worker_count, mp_context=mp.get_context("spawn"))
        self._mtf_executor_workers = worker_count
        return self._mtf_executor

    def _effective_mtf_worker_count(self, ray_total: int) -> int:
        if self.mtf_worker_count <= 1:
            return 1
        cpu_total = max(1, int(os.cpu_count() or 1))
        if ray_total < 1024 or cpu_total <= 1:
            return 1
        return max(1, min(self.mtf_worker_count, cpu_total, ray_total // 1024))

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

        if isinstance(operand, MagnificationOperand):
            _, _, _, a, _, _, _, _, _, _, _, _, _ = system.Parax(operand.wavelength)
            value = float(a)
            residual = value - operand.target
            return OperandResult(
                name=operand.resolved_name(),
                value=value,
                target=operand.target,
                residual=residual,
                weighted=operand.weight * residual * residual,
                metadata={},
            )

        if isinstance(operand, EntrancePupilPositionOperand):
            value, metadata = self._pupil_position(system, operand, which="entrance")
            residual = value - operand.target
            return OperandResult(
                name=operand.resolved_name(),
                value=value,
                target=operand.target,
                residual=residual,
                weighted=operand.weight * residual * residual,
                metadata=metadata,
            )

        if isinstance(operand, ExitPupilPositionOperand):
            value, metadata = self._pupil_position(system, operand, which="exit")
            residual = value - operand.target
            return OperandResult(
                name=operand.resolved_name(),
                value=value,
                target=operand.target,
                residual=residual,
                weighted=operand.weight * residual * residual,
                metadata=metadata,
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

        if isinstance(operand, MTFAtFrequencyOperand):
            value, metadata = self._mtf_at_frequency(system, operand)
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
        pupil = Kos.PupilCalc(
            system,
            operand.surface_index,
            operand.wavelength,
            "EPD",
            max(
                0.01,
                2.0
                * max((float(surface.Diameter) / 2.0 for surface in system.SDT), default=1.0)
                * max(min(operand.ray_height_factor, 1.5), 0.05),
            ),
        )
        pupil.Samp = max(2, int(operand.ray_count))
        pupil.Ptype = "hexapolar"
        pupil.FieldType = operand.field_type
        pupil.FieldX = float(operand.field_x)
        pupil.FieldY = float(operand.field_y)
        x, y, z, L, M, N = pupil.Pattern2Field()
        Kos.TraceLoop(x, y, z, L, M, N, operand.wavelength, rays, clean=1)
        X, Y, Z, L, M, N = rays.pick(operand.surface_index)
        if X.size == 0:
            raise RuntimeError("No valid rays reached the requested surface")
        rms, cen_x, cen_y = Kos.RMS(X, Y, Z, L, M, N)
        return float(rms), {
            "centroid_x": float(cen_x),
            "centroid_y": float(cen_y),
            "surface_index": operand.surface_index,
            "ray_count": int(pupil.Samp),
            "field_type": operand.field_type,
            "field_y": float(operand.field_y),
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

    def _pupil_position(self, system, operand, which: str):
        pupil = Kos.PupilCalc(
            system,
            operand.surface_index,
            operand.wavelength,
            operand.aperture_type,
            operand.aperture_value,
        )
        if which == "entrance":
            value = float(pupil.PosPupInp[2])
            radius = float(pupil.RadPupInp)
        else:
            value = float(pupil.PosPupOut[2])
            radius = float(pupil.RadPupOut)
        return value, {
            "surface_index": operand.surface_index,
            "aperture_type": operand.aperture_type,
            "aperture_value": float(operand.aperture_value),
            "radius": radius,
            "which": which,
        }

    def _mtf_at_frequency(self, system, operand: MTFAtFrequencyOperand):
        surface_index = max(0, int(operand.surface_index))
        aperture_type = str(getattr(operand, "aperture_type", "EPD")).strip().upper()
        if aperture_type not in {"STOP", "EPD"}:
            aperture_type = "EPD"
        aperture_value = float(getattr(operand, "aperture_value", 0.0))
        if aperture_value <= 0.0:
            aperture_value = max(
                0.01,
                2.0 * max((float(surface.Diameter) / 2.0 for surface in system.SDT), default=1.0),
            )

        rays = Kos.raykeeper(system)
        pupil = Kos.PupilCalc(
            system,
            surface_index,
            operand.wavelength,
            aperture_type,
            aperture_value,
        )
        pupil.Samp = max(2, int(operand.ray_count))
        pupil.Ptype = "hexapolar"
        pupil.FieldType = operand.field_type
        pupil.FieldX = float(operand.field_x)
        pupil.FieldY = float(operand.field_y)
        x, y, z, L, M, N = pupil.Pattern2Field()
        ray_total = int(np.asarray(x).size)
        worker_count = self._effective_mtf_worker_count(ray_total)
        x_local_parts: list[np.ndarray] = []
        y_local_parts: list[np.ndarray] = []
        if worker_count > 1:
            surface_specs = [self._surface_to_dict(surface) for surface in system.SDT]
            executor = self._ensure_mtf_executor(worker_count)
            if executor is not None:
                try:
                    futures = []
                    indices = np.array_split(np.arange(ray_total), worker_count)
                    for chunk in indices:
                        if chunk.size == 0:
                            continue
                        futures.append(
                            executor.submit(
                                _trace_mtf_chunk,
                                surface_specs,
                                float(operand.wavelength),
                                np.asarray(x, dtype=float)[chunk],
                                np.asarray(y, dtype=float)[chunk],
                                np.asarray(z, dtype=float)[chunk],
                                np.asarray(L, dtype=float)[chunk],
                                np.asarray(M, dtype=float)[chunk],
                                np.asarray(N, dtype=float)[chunk],
                            )
                        )
                    for future in futures:
                        x_part, y_part = future.result()
                        if x_part.size:
                            x_local_parts.append(x_part)
                            y_local_parts.append(y_part)
                except Exception:
                    x_local_parts.clear()
                    y_local_parts.clear()
                    worker_count = 1

        if worker_count <= 1:
            Kos.TraceLoop(x, y, z, L, M, N, operand.wavelength, rays, clean=1)
            X, Y, _Z, _L, _M, _N = _pick_image_plane_data_static(rays)
            x_local = np.asarray(X, dtype=float)
            y_local = np.asarray(Y, dtype=float)
        else:
            if not x_local_parts:
                raise RuntimeError("MTF parallel trace returned no image-plane samples")
            x_local = np.concatenate(x_local_parts)
            y_local = np.concatenate(y_local_parts)
        finite = np.isfinite(x_local) & np.isfinite(y_local)
        x_local = x_local[finite]
        y_local = y_local[finite]
        if x_local.size < 4:
            raise RuntimeError("Not enough image-plane ray samples for MTF")
        span_x = max(float(np.ptp(x_local)), 1e-3)
        span_y = max(float(np.ptp(y_local)), 1e-3)
        span = max(span_x, span_y) * 1.25
        bins = 128
        hist, xedges, yedges = np.histogram2d(
            x_local,
            y_local,
            bins=bins,
            range=[[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]],
        )
        psf = hist / max(np.sum(hist), 1.0)
        otf = np.fft.fftshift(np.fft.fft2(psf))
        mtf = np.abs(otf)
        mtf /= max(float(np.max(mtf)), 1e-12)
        dx = float(xedges[1] - xedges[0])
        freq = np.fft.fftshift(np.fft.fftfreq(bins, d=dx))
        center = bins // 2
        positive = freq[center:]
        tangential = np.asarray(mtf[center, center:], dtype=float)
        sagittal = np.asarray(mtf[center:, center], dtype=float)
        count = min(len(positive), len(tangential), len(sagittal))
        positive = positive[:count]
        tangential = tangential[:count]
        sagittal = sagittal[:count]
        mtf_mode = str(getattr(operand, "mode", "average")).strip().lower()
        if mtf_mode == "tangential":
            curve = tangential
        elif mtf_mode == "sagittal":
            curve = sagittal
        else:
            curve = 0.5 * (tangential + sagittal)
            mtf_mode = "average"
        value = float(np.interp(float(operand.frequency), positive, curve, left=curve[0], right=curve[-1]))
        return value, {
            "frequency": float(operand.frequency),
            "mode": mtf_mode,
            "ray_count": int(x_local.size),
            "worker_count": int(worker_count),
            "surface_index": surface_index,
            "aperture_type": aperture_type,
            "aperture_value": float(aperture_value),
            "field_type": operand.field_type,
            "field_y": float(operand.field_y),
        }

    @staticmethod
    def _surface_to_dict(surface) -> dict:
        return {
            "Name": getattr(surface, "Name", ""),
            "name": getattr(surface, "Name", ""),
            "Rc": float(getattr(surface, "Rc", 0.0)),
            "rc": float(getattr(surface, "Rc", 0.0)),
            "Thickness": float(getattr(surface, "Thickness", 0.0)),
            "thickness": float(getattr(surface, "Thickness", 0.0)),
            "Diameter": float(getattr(surface, "Diameter", 0.0)),
            "diameter": float(getattr(surface, "Diameter", 0.0)),
            "Glass": getattr(surface, "Glass", "AIR"),
            "glass": getattr(surface, "Glass", "AIR"),
            "Thin_Lens": float(getattr(surface, "Thin_Lens", 0.0)),
            "thin_lens": float(getattr(surface, "Thin_Lens", 0.0)),
            "diff_ord": float(getattr(surface, "Diff_Ord", 0.0)),
            "grating_d": float(getattr(surface, "Grating_D", 0.0)),
            "tilt_x": float(getattr(surface, "TiltX", 0.0)),
            "tilt_y": float(getattr(surface, "TiltY", 0.0)),
            "tilt_z": float(getattr(surface, "TiltZ", 0.0)),
            "desp_x": float(getattr(surface, "DespX", 0.0)),
            "desp_y": float(getattr(surface, "DespY", 0.0)),
            "desp_z": float(getattr(surface, "DespZ", 0.0)),
            "axis_move": float(getattr(surface, "AxisMove", 0.0)),
            "drawing": float(getattr(surface, "Drawing", 1.0)),
        }
