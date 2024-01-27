from __future__ import annotations

import numpy as np
import cmath

from abc import ABC
from typing import Callable
from logzero import logger
from dataclasses import dataclass

from sampling_fo2.fol.syntax import Pred
from sampling_fo2.utils import Rational
from sampling_fo2.utils.polynomial import coeff_dict, create_vars, Symbol, expand
from sampling_fo2.utils.third_typing import RingElement


class Constraint(ABC):
    pass


@dataclass(frozen=True)
class TreeConstraint(Constraint):
    pred: Pred

    def __str__(self):
        return "Tree({})".format(self.pred)

    def __repr__(self):
        return str(self)


class CardinalityConstraint(Constraint):
    def __init__(self, pred2card: dict[Pred, tuple[str, int]] = None):
        self.pred2card: dict[Pred, tuple[str, int]] = dict()
        if pred2card:
            self.pred2card: dict[Pred, tuple[str, int]] = pred2card


    def preds(self):
        return list(self.pred2card.keys())

    def transform_weighting(self, get_weight: Callable[[Pred], tuple[Rational, Rational]]) \
            -> dict[Pred, tuple[Rational, Rational]]:
        new_weights: dict[Pred, tuple[RingElement, RingElement]] = {}
        self.gen_vars = create_vars('x0:{}'.format(
            len(self.pred2card))
        )
        for sym, pred in zip(self.gen_vars, self.pred2card.keys()):
            weight = get_weight(pred)
            new_weights[pred] = (weight[0] * sym, weight[1])
        return new_weights

    def decode_poly(self, poly: RingElement) -> RingElement:
        poly = expand(poly)
        coeffs = coeff_dict(poly, self.gen_vars)
        # logger.debug('coeffs: %s', list(coeffs))
        res = Rational(0, 1)
        for degrees, coeff in coeffs:
            if self.valid(degrees):
                res += coeff
        return res

    def valid(self, degrees: list[int]) -> bool:
        return eval(self.validator.format(*degrees))

    def add(self, k: Pred, v: tuple[str, int]):
        self.pred2card[k] = v

    def build(self):
        self.gen_vars: list[Symbol] = list()
        self.validator: str = ''
        validator_list: list[str] = []
        for _, (comp, card) in self.pred2card.items():
            if comp == '=':
                comp = '=='
            validator_list.append('{{}} {} {}'.format(comp, card))
        self.validator = ' and '.join(validator_list)
        logger.info('cardinality validator: %s', self.validator)

    def __str__(self):
        s = ''
        for pred, (op, card) in self.pred2card.items():
            s += '|{}| {} {} and '.format(pred, op, card)
        return s

    def __repr__(self):
        return str(self)
