"""Optimization utilities for KrakenOS."""

from .specs import OPERAND_REGISTRY, VARIABLE_REGISTRY, OperandSpec, VariableSpec
from .variables import OpticalVariable
from .operands import (
    MeritOperand,
    SpotRMSOperand,
    EffectiveFocalLengthOperand,
    MagnificationOperand,
    MTFAtFrequencyOperand,
    EntrancePupilPositionOperand,
    ExitPupilPositionOperand,
    WavefrontRMSOperand,
    ThicknessPenaltyOperand,
    InvalidTracePenaltyOperand,
)
from .merit import MeritFunction, MeritResult, OperandResult
from .evaluator import MeritEvaluator
