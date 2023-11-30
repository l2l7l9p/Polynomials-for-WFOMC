from __future__ import annotations

import pandas as pd
import functools

from typing import FrozenSet, List, Tuple
from dataclasses import dataclass, field
from logzero import logger
from sympy import Poly

from sampling_fo2.fol.syntax import AtomicFormula, Pred, Term, a, b, X
from sampling_fo2.fol.utils import get_predicates
from sampling_fo2.utils import Rational


@dataclass(frozen=True)
class Cell(object):
    """
    In other words, the Unary types
    """
    code: Tuple[bool] = field(hash=False, compare=False)
    preds: Tuple[Pred] = field(hash=False, compare=False)
    # for hashing
    _identifier: FrozenSet[Tuple[Pred, bool]] = field(
        default=None, repr=False, init=False, hash=True, compare=True)

    def __post_init__(self):
        object.__setattr__(self, '_identifier',
                           frozenset(zip(self.preds, self.code)))

    @functools.lru_cache(maxsize=None)
    def get_evidences(self, term: Term) -> FrozenSet[AtomicFormula]:
        evidences: set[AtomicFormula] = set()
        for i, p in enumerate(self.preds):
            atom = p(*([term] * p.arity))
            if (self.code[i]):
                evidences.add(atom)
            else:
                evidences.add(~atom)
        return frozenset(evidences)

    @functools.lru_cache(maxsize=None)
    def is_positive(self, pred: Pred) -> bool:
        return self.code[self.preds.index(pred)]

    def negate(self, pred: Pred) -> Cell:
        idx = self.preds.index(pred)
        new_code = list(self.code)
        new_code[idx] = not new_code[idx]
        return Cell(tuple(new_code), self.preds)

    def drop_preds(self, preds: List[Pred] = None, prefixes: List[str] = None) -> Cell:
        if not preds and not prefixes:
            raise RuntimeError(
                'Dropped pred is not assigned'
            )
        if preds is not None:
            all_preds = [pred for pred in preds]
        else:
            all_preds = []

        if prefixes is not None:
            for prefix in prefixes:
                all_preds.extend(
                    get_predicates(prefix)
                )
        new_code, new_preds = zip(
            *[(c, p) for c, p in zip(self.code, self.preds)
              if p not in all_preds]
        )
        return Cell(tuple(new_code), tuple(new_preds))

    def __str__(self):
        evidences: frozenset[AtomicFormula] = self.get_evidences(X)
        # evidences = filter(lambda x: x.pred.name.startswith('skolem') or x.pred.name.startswith('aux') or x.pred.name == 'E', evidences)
        lits = [str(lit) for lit in evidences]
        lits.sort()
        return '^'.join(lits)

    def __repr__(self):
        return self.__str__()


class TwoTable(object):
    def __init__(self, model_table: pd.DataFrame, cell_1: Cell, cell_2: Cell):
        """
        The table containing the weight of all B-types.
        Note the order of cell_1 and cell_2 matters!
        """
        self.model_table = model_table
        self.cell_1: Cell = cell_1
        self.cell_2: Cell = cell_2

        self.evidences: FrozenSet[AtomicFormula] = frozenset(
            self.cell_1.get_evidences(a).union(
                self.cell_2.get_evidences(b)
            )
        )
        self.model_table = self._conditional_on(self.evidences)

    def _conditional_on(self, evidences: FrozenSet[AtomicFormula] = None) -> pd.DataFrame:
        if evidences is None:
            return self.model_table
        table = self.model_table
        for e in evidences:
            atom = e.make_positive()
            if atom in table.columns:
                if e.positive:
                    table = table[table[atom]]
                else:
                    table = table[~table[atom]]
                if table.empty:
                    return table
        return table

    def get_weight(self, evidences: FrozenSet[AtomicFormula] = None) -> Poly:
        if not self.satisfiable(evidences):
            return Rational(0, 1)
        table = self._conditional_on(evidences)
        return functools.reduce(
            lambda a, b: a + b,
            table.weight
        )

    def get_two_tables(self, evidences: FrozenSet[AtomicFormula] = None) \
            -> Tuple[FrozenSet[AtomicFormula], Poly]:
        two_tables = []
        df = self._conditional_on(evidences)
        if len(df) == 0:
            logger.warning(
                'Cell pair (%s, %s) with evidences %s is not satisfiable',
                self.cell_1, self.cell_2, evidences
            )
            return two_table
        for r in df.iterrows():
            two_table = set()
            for k, v in r[1].items():
                if k == 'weight':
                    weight = v
                else:
                    if v:
                        two_table.add(k)
                    else:
                        two_table.add(~k)
            two_tables.append((frozenset(two_table), weight))
        return two_tables

    def satisfiable(self, evidences: FrozenSet[AtomicFormula] = None) -> bool:
        table = self._conditional_on(evidences)
        if table.empty:
            return False
        return True
