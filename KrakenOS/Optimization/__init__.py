"""Optimization utilities for KrakenOS."""

from .variables import OpticalVariable
from .operands import (
    MeritOperand,
    SpotRMSOperand,
    EffectiveFocalLengthOperand,
    WavefrontRMSOperand,
    ThicknessPenaltyOperand,
    InvalidTracePenaltyOperand,
)
from .merit import MeritFunction, MeritResult, OperandResult
from .evaluator import MeritEvaluator
