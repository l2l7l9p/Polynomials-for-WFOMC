from __future__ import annotations

from sampling_fo2.fol.sc2 import SC2, to_sc2
from sampling_fo2.fol.syntax import Const, Pred, top, AUXILIARY_PRED_NAME, \
    Formula, QuantifiedFormula, Universal, Equivalence
from sampling_fo2.fol.utils import new_predicate
from sampling_fo2.network.constraint import CardinalityConstraint
from sampling_fo2.utils.polynomial import Rational
from fractions import Fraction
import math


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


class MLNProblem(object):
    """
    A Markov Logic Network problem.
    """

    def __init__(self, rules: tuple[list[tuple[Rational, Rational]], list[Formula]],
                 domain: set[Const],
                 cardinality_constraint: CardinalityConstraint):
        self.rules = rules
        # self.formulas: rules[1]
        # self.formula_weights: = dict(zip(rules[1], rules[0]))
        self.domain: set[Const] = domain
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint


def MLN_to_WFOMC(mln: MLNProblem):
    sentence = top
    weightings: dict[Pred, tuple[Rational, Rational]] = dict()
    for weighting, formula in zip(*mln.rules):
        free_vars = formula.free_vars()
        if weighting != float('inf'):
            aux_pred = new_predicate(len(free_vars), AUXILIARY_PRED_NAME)
            formula = Equivalence(formula, aux_pred(*free_vars))
            weightings[aux_pred] = (Rational(Fraction(math.exp(weighting)).numerator,
                                             Fraction(math.exp(weighting)).denominator), Rational(1, 1))
        for free_var in free_vars:
            formula = QuantifiedFormula(Universal(free_var), formula)
        sentence = sentence & formula

    sc2 = to_sc2(sentence)
    return WFOMCSProblem(sc2, mln.domain, weightings, mln.cardinality_constraint)
