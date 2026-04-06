from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

from .merit import MeritFunction
from .operands import (
    EffectiveFocalLengthOperand,
    EntrancePupilPositionOperand,
    ExitPupilPositionOperand,
    MagnificationOperand,
    MTFAtFrequencyOperand,
    SpotRMSOperand,
    ThicknessPenaltyOperand,
    WavefrontRMSOperand,
)


@dataclass(frozen=True)
class VariableSpec:
    """Metadata describing one editable optimization variable type."""

    key: str
    label: str
    field: str
    parameter: str
    optimize_attr: str
    bounds_attr: str
    supported_surfaces: tuple[str, ...]
    default_bounds: Callable[[float], tuple[float, float]]

    def is_supported(self, row) -> bool:
        return row.surface in self.supported_surfaces

    def is_enabled(self, row) -> bool:
        return bool(getattr(row, self.optimize_attr))

    def set_enabled(self, row, enabled: bool) -> None:
        setattr(row, self.optimize_attr, enabled)

    def get_bounds(self, row) -> tuple[float, float] | None:
        return getattr(row, self.bounds_attr)

    def set_bounds(self, row, bounds: tuple[float, float] | None) -> None:
        setattr(row, self.bounds_attr, bounds)

    def value_from_row(self, row) -> float:
        return float(getattr(row, self.field))


@dataclass(frozen=True)
class OperandSpec:
    """Metadata describing one merit-function option shown in the UI."""

    key: str
    label: str
    description: str
    default_weight: float
    default_target: float
    build_operands: Callable[[object], Iterable[object]]
    controls: tuple[str, ...] = ("weight", "target")

    def build_merit_function(self, editor) -> MeritFunction:
        return MeritFunction(list(self.build_operands(editor)))


def _default_rc_bounds(value: float) -> tuple[float, float]:
    if abs(value) < 1e-6:
        return (-100.0, 100.0)
    scale = max(abs(value) * 0.5, 5.0)
    return (value - scale, value + scale)


def _default_thickness_bounds(value: float) -> tuple[float, float]:
    if value <= 0.0:
        return (0.01, 10.0)
    lower = max(0.01, value * 0.5)
    upper = max(lower + 0.5, value * 1.5)
    return (lower, upper)


VARIABLE_REGISTRY: dict[str, VariableSpec] = {
    "rc": VariableSpec(
        key="rc",
        label="Radius",
        field="rc",
        parameter="Rc",
        optimize_attr="optimize_rc",
        bounds_attr="optimize_rc_bounds",
        supported_surfaces=("Standard",),
        default_bounds=_default_rc_bounds,
    ),
    "thickness": VariableSpec(
        key="thickness",
        label="Thickness",
        field="thickness",
        parameter="Thickness",
        optimize_attr="optimize_thickness",
        bounds_attr="optimize_thickness_bounds",
        supported_surfaces=("Standard", "Thin Lens", "Grating"),
        default_bounds=_default_thickness_bounds,
    ),
}


def _spot_rms_operands(editor) -> list[object]:
    return [
        SpotRMSOperand(
            name="Spot RMS",
            weight=editor._operand_weight("Spot RMS"),
            target=editor._operand_target("Spot RMS"),
            surface_index=-1,
            wavelength=editor._operand_wavelength("Spot RMS"),
            ray_count=max(5, editor._current_ray_count()),
            ray_height_factor=editor._current_ray_height_factor(),
            field_type=editor._operand_field_type("Spot RMS"),
            field_y=editor._operand_field("Spot RMS"),
        )
    ]


def _wavefront_rms_operands(editor) -> list[object]:
    return [
        WavefrontRMSOperand(
            name="Wavefront RMS",
            weight=editor._operand_weight("Wavefront RMS"),
            target=editor._operand_target("Wavefront RMS"),
            surface_index=editor._analysis_surface_index(),
            wavelength=editor._operand_wavelength("Wavefront RMS"),
            aperture_type=editor._current_aperture_type(),
            aperture_value=editor._current_aperture_value(),
            sample_size=9,
            field_type=editor._operand_field_type("Wavefront RMS"),
            field_y=editor._operand_field("Wavefront RMS"),
        )
    ]


def _effective_focal_length_operands(editor) -> list[object]:
    return [
        EffectiveFocalLengthOperand(
            name="EFFL",
            weight=editor._operand_weight("EFFL"),
            target=editor._operand_target("EFFL"),
            wavelength=editor._operand_wavelength("EFFL"),
        )
    ]


def _magnification_operands(editor) -> list[object]:
    return [
        MagnificationOperand(
            name="Magnification",
            weight=editor._operand_weight("Magnification"),
            target=editor._operand_target("Magnification"),
            wavelength=editor._operand_wavelength("Magnification"),
        )
    ]


def _entrance_pupil_position_operands(editor) -> list[object]:
    return [
        EntrancePupilPositionOperand(
            name="Entrance pupil z",
            weight=editor._operand_weight("Entrance pupil z"),
            target=editor._operand_target("Entrance pupil z"),
            surface_index=editor._analysis_surface_index(),
            wavelength=editor._operand_wavelength("Entrance pupil z"),
            aperture_type=editor._current_aperture_type(),
            aperture_value=editor._current_aperture_value(),
        )
    ]


def _exit_pupil_position_operands(editor) -> list[object]:
    return [
        ExitPupilPositionOperand(
            name="Exit pupil z",
            weight=editor._operand_weight("Exit pupil z"),
            target=editor._operand_target("Exit pupil z"),
            surface_index=editor._analysis_surface_index(),
            wavelength=editor._operand_wavelength("Exit pupil z"),
            aperture_type=editor._current_aperture_type(),
            aperture_value=editor._current_aperture_value(),
        )
    ]


def _thickness_penalty_operands(editor) -> list[object]:
    operands: list[object] = []
    minimum = editor._operand_target("Thickness penalty")
    weight = editor._operand_weight("Thickness penalty")
    for variable in editor._build_optimization_variables():
        if variable.parameter != "Thickness":
            continue
        operands.append(
            ThicknessPenaltyOperand(
                name=f"{variable.normalized_name()} >= {minimum:g}",
                weight=weight,
                target=0.0,
                surface_index=variable.surface_index,
                minimum=minimum,
            )
        )
    return operands


def _mtf_at_frequency_operands(editor) -> list[object]:
    return [
        MTFAtFrequencyOperand(
            name="MTF @ freq",
            weight=editor._operand_weight("MTF @ freq"),
            target=editor._operand_target("MTF @ freq"),
            surface_index=editor._operand_surface_index("MTF @ freq"),
            wavelength=editor._operand_wavelength("MTF @ freq"),
            frequency=editor._current_mtf_frequency(),
            ray_count=max(24, editor._current_ray_count() * 6),
            mode=editor._operand_mtf_mode("MTF @ freq"),
            field_type=editor._operand_field_type("MTF @ freq"),
            field_x=editor._operand_field_x("MTF @ freq"),
            field_y=editor._operand_field_y("MTF @ freq"),
        )
    ]

OPERAND_REGISTRY: dict[str, OperandSpec] = {
    "spot_rms": OperandSpec(
        key="spot_rms",
        label="Spot RMS",
        description="Minimize the RMS spot size on the image surface.",
        default_weight=1.0,
        default_target=0.0,
        build_operands=_spot_rms_operands,
        controls=("weight", "target", "wavelength", "field"),
    ),
    "wavefront_rms": OperandSpec(
        key="wavefront_rms",
        label="Wavefront RMS",
        description="Minimize sampled wavefront RMS on the selected analysis surface.",
        default_weight=1.0,
        default_target=0.0,
        build_operands=_wavefront_rms_operands,
        controls=("weight", "target", "wavelength", "field", "surface", "aperture", "aperture_value"),
    ),
    "effective_focal_length": OperandSpec(
        key="effective_focal_length",
        label="EFFL",
        description="Drive the system effective focal length toward a target value.",
        default_weight=1.0,
        default_target=100.0,
        build_operands=_effective_focal_length_operands,
        controls=("weight", "target", "wavelength"),
    ),
    "magnification": OperandSpec(
        key="magnification",
        label="Magnification",
        description="Drive the paraxial transverse magnification toward a target value.",
        default_weight=1.0,
        default_target=1.0,
        build_operands=_magnification_operands,
        controls=("weight", "target", "wavelength"),
    ),
    "entrance_pupil_position": OperandSpec(
        key="entrance_pupil_position",
        label="Entrance pupil z",
        description="Drive the entrance pupil axial position toward a target value.",
        default_weight=1.0,
        default_target=0.0,
        build_operands=_entrance_pupil_position_operands,
        controls=("weight", "target", "wavelength", "surface", "aperture", "aperture_value"),
    ),
    "exit_pupil_position": OperandSpec(
        key="exit_pupil_position",
        label="Exit pupil z",
        description="Drive the exit pupil axial position toward a target value.",
        default_weight=1.0,
        default_target=0.0,
        build_operands=_exit_pupil_position_operands,
        controls=("weight", "target", "wavelength", "surface", "aperture", "aperture_value"),
    ),
    "thickness_penalty": OperandSpec(
        key="thickness_penalty",
        label="Thickness penalty",
        description="Penalize optimized thickness variables that fall below the target minimum.",
        default_weight=1.0,
        default_target=0.1,
        build_operands=_thickness_penalty_operands,
        controls=("weight", "target"),
    ),
    "mtf_at_frequency": OperandSpec(
        key="mtf_at_frequency",
        label="MTF @ freq",
        description="Drive the geometric MTF toward a target value at the selected spatial frequency.",
        default_weight=1.0,
        default_target=0.5,
        build_operands=_mtf_at_frequency_operands,
        controls=("weight", "target", "wavelength", "field_xy", "surface", "aperture", "aperture_value", "frequency", "mtf_mode"),
    ),
}
