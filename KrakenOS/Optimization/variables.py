from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class OpticalVariable:
    """A scalar design variable attached to a Kraken surface parameter."""

    surface_index: int
    parameter: str
    lower_bound: float
    upper_bound: float
    scale: float = 1.0
    enabled: bool = True
    name: str = ""

    def normalized_name(self) -> str:
        if self.name:
            return self.name
        return f"surf[{self.surface_index}].{self.parameter}"

