from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class OperandResult:
    name: str
    value: float
    target: float
    residual: float
    weighted: float
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class MeritResult:
    total: float
    operands: list[OperandResult]
    valid: bool
    message: str = ""
    prescription: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class MeritFunction:
    operands: list[Any]
    invalid_penalty: float = 1e9

