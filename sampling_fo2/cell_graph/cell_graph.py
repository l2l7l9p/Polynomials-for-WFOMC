from __future__ import annotations

import pandas as pd
import functools
import networkx as nx

from typing import Callable, Dict, FrozenSet, List, Tuple
from logzero import logger
from sympy import Poly
from copy import deepcopy

from sampling_fo2.fol.sc2 import SC2

from sampling_fo2.fol.syntax import AtomicFormula, Const, Pred, QFFormula, a, b, c
from sampling_fo2.utils import Rational, RingElement
from sampling_fo2.utils.multinomial import MultinomialCoefficients

from .components import Cell, TwoTable


class CellGraph(object):
    """
    Cell graph that handles cells (i.e., 1-types, in the sampling paper) and the wmc between them.
    """

    def __init__(self, formula: QFFormula,
                 get_weight: Callable[[Pred], Tuple[RingElement, RingElement]]):
        """
        Cell graph that handles cells (1-types) and the WMC between them

        :param sentence CNF: the sentence in the form of CNF
        :param get_weight Callable[[Pred], Tuple[mpq, mpq]]: the weighting function
        :param conditional_formulas List[CNF]: the optional conditional formula appended in WMC computing
        """
        self.formula: QFFormula = formula
        self.get_weight: Callable[[Pred],
                                  Tuple[RingElement, RingElement]] = get_weight
        self.preds: Tuple[Pred] = tuple(self.formula.preds())
        logger.debug('prednames: %s', self.preds)

        gnd_formula_ab1: QFFormula = self._ground_on_tuple(
            self.formula, a, b
        )
        gnd_formula_ab2: QFFormula = self._ground_on_tuple(
            self.formula, b, a
        )
        self.gnd_formula_ab: QFFormula = \
            gnd_formula_ab1 & gnd_formula_ab2
        self.gnd_formula_cc: QFFormula = self._ground_on_tuple(
            self.formula, c
        )
        logger.info('ground a b: %s', self.gnd_formula_ab)
        logger.info('ground c: %s', self.gnd_formula_cc)

        # build cells
        self.cells: List[Cell] = self._build_cells()
        # filter cells
        logger.info('the number of valid cells: %s',
                    len(self.cells))

        self.cell_weights: Dict[Cell, Poly] = self._compute_cell_weights()
        self.two_tables: Dict[Tuple[Cell, Cell],
                              TwoTable] = self._build_two_tables()

    def _ground_on_tuple(self, formula: QFFormula,
                         c1: Const, c2: Const = None) -> QFFormula:
        variables = formula.vars()
        if len(variables) > 2:
            raise RuntimeError(
                "Can only ground out FO2"
            )
        if len(variables) == 1:
            constants = [c1]
        else:
            if c2 is not None:
                constants = [c1, c2]
            else:
                constants = [c1, c1]
        substitution = dict(zip(variables, constants))
        gnd_formula = formula.substitute(substitution)
        return gnd_formula

    def show(self):
        logger.info(str(self))

    def __str__(self):
        s = 'CellGraph:\n'
        s += 'predicates: {}\n'.format(self.preds)
        cell_weight_df = []
        twotable_weight_df = []
        for _, cell1 in enumerate(self.cells):
            cell_weight_df.append(
                [str(cell1), self.get_cell_weight(cell1)]
            )
            twotable_weight = []
            for _, cell2 in enumerate(self.cells):
                # if idx1 < idx2:
                #     twotable_weight.append(0)
                #     continue
                twotable_weight.append(
                    self.get_two_table_weight(
                        (cell1, cell2))
                )
            twotable_weight_df.append(twotable_weight)
        cell_str = [str(cell) for cell in self.cells]
        cell_weight_df = pd.DataFrame(cell_weight_df, index=None,
                                      columns=['Cell', 'Weight'])
        twotable_weight_df = pd.DataFrame(twotable_weight_df, index=cell_str,
                                          columns=cell_str)
        s += 'cell weights: \n'
        s += cell_weight_df.to_markdown() + '\n'
        s += '2table weights: \n'
        s += twotable_weight_df.to_markdown()
        return s

    def __repr__(self):
        return str(self)

    def get_cells(self, cell_filter: Callable[[Cell], bool] = None) -> List[Cell]:
        if cell_filter is None:
            return self.cells
        return list(filter(cell_filter, self.cells))

    @functools.lru_cache(maxsize=None, typed=True)
    def get_cell_weight(self, cell: Cell) -> Poly:
        if cell not in self.cell_weights:
            logger.warning(
                "Cell %s not found", cell
            )
            return 0
        return self.cell_weights.get(cell)

    def _check_existence(self, cells: Tuple[Cell, Cell]):
        if cells not in self.two_tables:
            raise ValueError(
                "Cells (%s) not found, note that the order of cells matters!", cells
            )

    @functools.lru_cache(maxsize=None, typed=True)
    def get_two_table_weight(self, cells: Tuple[Cell, Cell],
                             evidences: FrozenSet[AtomicFormula] = None) -> RingElement:
        self._check_existence(cells)
        return self.two_tables.get(cells).get_weight(evidences)

    def get_all_weights(self) -> Tuple[List[RingElement], List[RingElement]]:
        cell_weights = []
        twotable_weights = []
        for i, cell_i in enumerate(self.cells):
            cell_weights.append(self.get_cell_weight(cell_i))
            twotable_weight = []
            for j, cell_j in enumerate(self.cells):
                if i > j:
                    twotable_weight.append(Rational(1, 1))
                else:
                    twotable_weight.append(self.get_two_table_weight(
                        (cell_i, cell_j)
                    ))
            twotable_weights.append(twotable_weight)
        return cell_weights, twotable_weights

    @functools.lru_cache(maxsize=None, typed=True)
    def satisfiable(self, cells: Tuple[Cell, Cell],
                    evidences: FrozenSet[AtomicFormula] = None) -> bool:
        self._check_existence(cells)
        return self.two_tables.get(cells).satisfiable(evidences)

    @functools.lru_cache(maxsize=None)
    def get_two_tables(self, cells: Tuple[Cell, Cell],
                       evidences: FrozenSet[AtomicFormula] = None) \
            -> Tuple[FrozenSet[AtomicFormula], RingElement]:
        self._check_existence(cells)
        return self.two_tables.get(cells).get_two_tables(evidences)

    def _build_cells(self):
        cells = []
        code = {}
        for model in self.gnd_formula_cc.models():
            for lit in model:
                code[lit.pred] = lit.positive
            cells.append(Cell(tuple(code[p] for p in self.preds), self.preds))
        return cells

    def _compute_cell_weights(self):
        weights = dict()
        for cell in self.cells:
            weight = Rational(1, 1)
            for i, pred in zip(cell.code, cell.preds):
                if pred.arity > 0:
                    if i:
                        weight = weight * self.get_weight(pred)[0]
                    else:
                        weight = weight * self.get_weight(pred)[1]
            weights[cell] = weight
        return weights

    @functools.lru_cache(maxsize=None)
    def get_nullary_weight(self, cell: Cell) -> RingElement:
        weight = Rational(1, 1)
        for i, pred in zip(cell.code, cell.preds):
            if pred.arity == 0:
                if i:
                    weight = weight * self.get_weight(pred)[0]
                else:
                    weight = weight * self.get_weight(pred)[1]
        return weight

    def _build_two_tables(self):
        # build a pd.DataFrame containing all model as well as the weight
        atoms = list(self.gnd_formula_ab.atoms())
        table = []
        for model in self.gnd_formula_ab.models():
            weight = Rational(1, 1)
            values = {}
            for lit in model:
                values[lit.make_positive()] = lit.positive
                # ignore the weight appearing in cell weight
                if (not (len(lit.args) == 1 or all(arg == lit.args[0]
                                                   for arg in lit.args))):
                    weight *= (self.get_weight(lit.pred)[0] if lit.positive else
                               self.get_weight(lit.pred)[1])
            table.append([values[atom] for atom in atoms] + [weight])
        model_table = pd.DataFrame(
            table, columns=atoms + ['weight']
        )
        # build twotable tables
        tables = dict()
        for cell in self.cells:
            for other_cell in self.cells:
                tables[(cell, other_cell)] = TwoTable(
                    model_table, cell, other_cell
                )
        return tables


class OptimizedCellGraph(CellGraph):
    def __init__(self, formula: QFFormula,
                 get_weight: Callable[[Pred], Tuple[RingElement, RingElement]],
                 domain_size: int):
        """
        Optimized cell graph for FastWFOMC
        :param formula: the formula to be grounded
        :param get_weight: a function that returns the weight of a predicate
        :param domain_size: the domain size
        """
        super().__init__(formula, get_weight)
        self.domain_size: int = domain_size
        self.cliques: list[list[Cell]] = self.build_symmetric_cliques()
        MultinomialCoefficients.setup(self.domain_size)
        self.i1_ind, self.i2_ind, \
            self.ind, self.nonind = self.find_independent_cliques()
        self.nonind_map: dict[int, int] = dict(
            zip(self.nonind, range(len(self.nonind))))

        self.term_cache = dict()

    def build_symmetric_cliques(self) -> List[List[Cell]]:
        cliques: list[list[Cell]] = []
        cells = deepcopy(self.get_cells())
        while len(cells) > 0:
            cell = cells.pop()
            clique = [cell]
            for other_cell in cells:
                if self._matches(clique, other_cell):
                    clique.append(other_cell)
            for other_cell in clique[1:]:
                cells.remove(other_cell)
            cliques.append(clique)
        cliques.sort(key=len)
        logger.info("Built %s symmetric cliques: %s", len(cliques), cliques)
        return cliques

    def find_independent_cliques(self) -> tuple[list[int], list[int], list[int], list[int]]:
        g = nx.Graph()
        g.add_nodes_from(range(len(self.cliques)))
        for i in range(len(self.cliques)):
            for j in range(i + 1, len(self.cliques)):
                if self.get_two_table_weight(
                        (self.cliques[i][0], self.cliques[j][0])
                ) != Rational(1, 1):
                    g.add_edge(i, j)

        self_loop = set()
        for i in range(len(self.cliques)):
            for j in range(self.domain_size):
                if self.get_J_term(i, j) != Rational(1, 1):
                    self_loop.add(i)
                    break

        g_ind = set(nx.maximal_independent_set(g))
        i2_ind = g_ind.intersection(self_loop)
        i1_ind = g_ind.difference(i2_ind)
        non_ind = g.nodes - i1_ind - i2_ind
        logger.info("Found i1 independent cliques: %s", i1_ind)
        logger.info("Found i2 independent cliques: %s", i2_ind)
        logger.info("Found non-independent cliques: %s", non_ind)
        return list(i1_ind), list(i2_ind), list(g_ind), list(non_ind)

    def _matches(self, clique, other_cell) -> bool:
        cell = clique[0]
        if self.get_cell_weight(cell) != self.get_cell_weight(other_cell) or \
                self.get_two_table_weight((cell, cell)) != self.get_two_table_weight((other_cell, other_cell)):
            return False

        if len(clique) > 1:
            third_cell = clique[1]
            r = self.get_two_table_weight((cell, third_cell))
            for third_cell in clique:
                if r != self.get_two_table_weight((other_cell, third_cell)):
                    return False

        for third_cell in self.get_cells():
            if other_cell == third_cell or third_cell in clique:
                continue
            r = self.get_two_table_weight((cell, third_cell))
            if r != self.get_two_table_weight((other_cell, third_cell)):
                return False
        return True

    def setup_term_cache(self):
        self.term_cache = dict()

    def get_term(self, iv: int, bign: int, partition: tuple[int]) -> RingElement:
        if (iv, bign) in self.term_cache:
            return self.term_cache[(iv, bign)]

        if iv == 0:
            accum = Rational(0, 1)
            for j in self.i1_ind:
                tmp = self.get_cell_weight(self.cliques[j][0])
                for i in self.nonind:
                    tmp = tmp * self.get_two_table_weight(
                        (self.cliques[i][0], self.cliques[j][0])) ** partition[self.nonind_map[i]]
                accum = accum + tmp
            accum = accum ** (self.domain_size - sum(partition) - bign)
            self.term_cache[(iv, bign)] = accum
            return accum
        else:
            sumtoadd = 0
            s = self.i2_ind[len(self.i2_ind) - iv]
            for nval in range(self.domain_size - sum(partition) - bign + 1):
                smul = MultinomialCoefficients.comb(
                    self.domain_size - sum(partition) - bign, nval
                )
                smul = smul * self.get_J_term(s, nval)
                smul = smul * self.get_cell_weight(self.cliques[s][0]) ** nval

                for i in self.nonind:
                    smul = smul * self.get_two_table_weight(
                        (self.cliques[i][0], self.cliques[s][0])
                    ) ** (partition[self.nonind_map[i]] * nval)
                smul = smul * self.get_term(
                    iv - 1, bign + nval, partition
                )
                sumtoadd = sumtoadd + smul
            self.term_cache[(iv, bign)] = sumtoadd
            return sumtoadd

    @functools.lru_cache(maxsize=None)
    def get_J_term(self, l: int, nhat: int) -> RingElement:
        if len(self.cliques[l]) == 1:
            thesum = self.get_two_table_weight(
                (self.cliques[l][0], self.cliques[l][0])
            ) ** (int(nhat * (nhat - 1) / 2))
        else:
            r = self.get_two_table_weight(
                (self.cliques[l][0], self.cliques[l][1]))
            thesum = (
                (r ** MultinomialCoefficients.comb(nhat, 2)) *
                self.get_d_term(l, 1, len(self.cliques[l]), nhat)
            )
        return thesum

    @functools.lru_cache(maxsize=None)
    def get_d_term(self, l: int, cur: int, maxi: int, n: int) -> RingElement:
        r = self.get_two_table_weight((self.cliques[l][0], self.cliques[l][1]))
        s = self.get_two_table_weight((self.cliques[l][0], self.cliques[l][0]))
        if cur == maxi:
            ret = (s / r) ** MultinomialCoefficients.comb(n, 2)
        else:
            ret = 0
            for ni in range(n + 1):
                mult = MultinomialCoefficients.comb(n, ni)
                mult = mult * ((s / r) ** MultinomialCoefficients.comb(ni, 2))
                mult = mult * self.get_d_term(l, cur + 1, maxi, n - ni)
                ret = ret + mult
        return ret
