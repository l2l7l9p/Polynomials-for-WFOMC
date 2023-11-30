from __future__ import annotations

from sampling_fo2.fol.sc2 import SC2
from sampling_fo2.fol.syntax import Const, Pred
from sampling_fo2.network.constraint import CardinalityConstraint
from sampling_fo2.utils.polynomial import Rational


class WFOMCSProblem(object):
    """
    A weighted first-order model counting/sampling problem.
    """

    def __init__(self, sentence: SC2,
                 domain: set[Const],
                 weights: dict[Pred, tuple[Rational, Rational]],
                 cardinality_constraint: CardinalityConstraint):
        self.domain: set[Const] = domain
        self.sentence: SC2 = sentence
        self.weights: dict[Pred, tuple[Rational, Rational]] = weights
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint
