import numpy as np
import cmath

from abc import ABC
from typing import FrozenSet, List, Tuple
from logzero import logger
from itertools import product
from copy import deepcopy
from dataclasses import dataclass

from sampling_fo2.fol.syntax import Pred, QuantifiedFormula, x, y
from sampling_fo2.network.mln import MLN, ComplexMLN


class Constraint(ABC):
    pass


@dataclass(frozen=True)
class TreeConstraint(Constraint):
    pred: Pred

    def __str__(self):
        return "Tree({})".format(self.pred)

    def __repr__(self):
        return str(self)


@dataclass(frozen=True)
class CardinalityConstraint(Constraint):
    pred2card: FrozenSet[Tuple[Pred, int]]

    def preds(self):
        return list(self.pred2card.keys())

    def valid(self, cards: List[int]):
        return cards == tuple(self.pred2card.values())

    def dft(self, mln: MLN) -> Tuple[ComplexMLN, np.ndarray, np.ndarray]:
        new_formulas = []
        dft_domain = []
        top_weights = []
        M = []
        for pred in self.preds():
            cnf = QuantifiedFormula.from_pred(pred, [x, y])
            new_formulas.append(cnf)
            D_f = mln.domain_size() ** pred.arity
            dft_domain.append(range(D_f + 1))
            M.append(D_f + 1)

        new_weights = [[] for _ in range(len(new_formulas))]
        M = np.array(M)
        logger.debug('dft domain: %s', dft_domain)
        d = np.prod(M)
        for i in product(*dft_domain):
            weight_for_constraint_formulas = complex(
                0, -2 * cmath.pi
            ) * np.array(i) / M
            for k, f in enumerate(new_formulas):
                new_weights[k].append(weight_for_constraint_formulas[k])
            if not self.valid(i):
                continue
            top_w = []
            for j in product(*dft_domain):
                top_w.append(cmath.exp(
                    complex(0, 2 * cmath.pi *
                            np.dot(np.array(i), np.array(j) / M))
                ))
            top_weights.append(top_w)
        top_weights = np.array(top_weights, dtype=np.complex256)
        new_weights = [np.array(w, dtype=np.complex256) for w in new_weights]
        formulas = mln.formulas + new_formulas
        weights = [
            np.tile(np.complex256(weight), int(d)) for weight in mln.weights
        ]
        weights.extend(new_weights)
        return ComplexMLN(formulas, weights, deepcopy(mln.domain)), top_weights, d

    def __str__(self):
        s = ''
        for pred, card in self.pred2card:
            s += '|{}| = {}\n'.format(pred, card)
        return s

    def __repr__(self):
        return str(self)
