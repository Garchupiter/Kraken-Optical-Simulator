from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MeritOperand:
    weight: float = 1.0
    target: float = 0.0
    enabled: bool = True
    name: str = ""

    def resolved_name(self) -> str:
        return self.name or self.__class__.__name__


@dataclass(frozen=True)
class SpotRMSOperand(MeritOperand):
    surface_index: int = -1
    wavelength: float = 0.55
    ray_count: int = 11
    ray_height_factor: float = 0.8
    field_type: str = "angle"
    field_x: float = 0.0
    field_y: float = 0.0


@dataclass(frozen=True)
class EffectiveFocalLengthOperand(MeritOperand):
    wavelength: float = 0.55


@dataclass(frozen=True)
class MagnificationOperand(MeritOperand):
    wavelength: float = 0.55


@dataclass(frozen=True)
class EntrancePupilPositionOperand(MeritOperand):
    surface_index: int = 1
    wavelength: float = 0.55
    aperture_type: str = "EPD"
    aperture_value: float = 10.0


@dataclass(frozen=True)
class ExitPupilPositionOperand(MeritOperand):
    surface_index: int = 1
    wavelength: float = 0.55
    aperture_type: str = "EPD"
    aperture_value: float = 10.0


@dataclass(frozen=True)
class WavefrontRMSOperand(MeritOperand):
    surface_index: int = 1
    wavelength: float = 0.55
    aperture_type: str = "EPD"
    aperture_value: float = 10.0
    sample_size: int = 9
    pattern: str = "hexapolar"
    field_type: str = "angle"
    field_x: float = 0.0
    field_y: float = 0.0


@dataclass(frozen=True)
class MTFAtFrequencyOperand(MeritOperand):
    surface_index: int = -1
    wavelength: float = 0.55
    frequency: float = 50.0
    ray_count: int = 30
    mode: str = "average"
    aperture_type: str = "EPD"
    aperture_value: float = 10.0
    field_type: str = "angle"
    field_x: float = 0.0
    field_y: float = 0.0


@dataclass(frozen=True)
class ThicknessPenaltyOperand(MeritOperand):
    surface_index: int = 0
    minimum: float = 0.0
    maximum: float | None = None


@dataclass(frozen=True)
class InvalidTracePenaltyOperand(MeritOperand):
    penalty: float = 1e9
