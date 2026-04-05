from __future__ import annotations

from dataclasses import dataclass

from ..evaluator import MeritEvaluator
from ..variables import OpticalVariable


@dataclass
class Pygmo2MeritProblem:
    """Thin adapter exposing a Kraken merit evaluator as a pygmo2 UDP."""

    evaluator: MeritEvaluator
    variables: list[OpticalVariable]

    def fitness(self, x):
        result = self.evaluator.evaluate(self.variables, x)
        return [result.total]

    def get_bounds(self):
        lower = [variable.lower_bound for variable in self.variables]
        upper = [variable.upper_bound for variable in self.variables]
        return (lower, upper)

    def get_name(self):
        return "Kraken Merit Function Problem"

    def get_extra_info(self):
        return "\n".join(variable.normalized_name() for variable in self.variables)

