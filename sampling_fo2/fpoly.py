'''
Input: an FO formula
Output: F-poly or Tutte poly of it
'''


from __future__ import annotations
from collections import defaultdict

from enum import Enum
import os
import argparse
import logging
import logzero

from logzero import logger
from typing import Callable
from contexttimer import Timer

from sampling_fo2.utils import PREDS_FOR_EXISTENTIAL, MultinomialCoefficients, multinomial, \
    multinomial_less_than, RingElement, Rational, round_rational, lagrange
from sampling_fo2.cell_graph import CellGraph, Cell
from sampling_fo2.context import WFOMCContext
from sampling_fo2.parser import parse_input
from sampling_fo2.fol.syntax import Const, Pred, QFFormula, AtomicFormula
from sampling_fo2.utils.polynomial import coeff_dict, create_vars, expand
from sampling_fo2.network.constraint import CardinalityConstraint


class Func(Enum):
    FPOLY = 'fpoly'
    TUTTE = 'tutte'

    def __str__(self):
        return self.value


def get_two_table_weight_with_v( cell_graph: CellGraph,
                                 cell_i: Cell,
                                 cell_j: Cell,
                                 special_pred : Pred,
                                 v: int = None ) :
    if (v==None) :
        return cell_graph.get_two_table_weight((cell_i, cell_j))
    else :
        evidence_special_pred_true = frozenset([
            AtomicFormula(special_pred,(Const('a'), Const('b')),True),
            AtomicFormula(special_pred,(Const('b'), Const('a')),True) ])
        weight_true = cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_true)
        evidence_special_pred_false = frozenset([
            AtomicFormula(special_pred,(Const('a'), Const('b')),False),
            AtomicFormula(special_pred,(Const('b'), Const('a')),False) ])
        weight_false = cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_false)
        return weight_true * v + weight_false


def F_poly_evaluate(cell_graph: CellGraph,
                    domain_size: int,
                    special_pred: Pred,
                    constraint: CardinalityConstraint,
                    v: int = None ) -> list :
    cells = cell_graph.get_cells()
    n_cells = len(cells)
    # for cell, n in cell_config.items():
        # if n > 0:
            # NOTE: nullary weight is multiplied once
            # val_cur = val_cur * cell_graph.get_nullary_weight(cell)
            # break
    dp_last : dict[tuple[int], Rational] = {tuple([0 for i in cells]) : 1}
    dp_cur : dict[tuple[int], Rational] = {}
    evaluate = []
    for u in range(domain_size+1) :
        for partition_last, val_last in dp_last.items() :
            domain_size_rest = domain_size
            for n_i in partition_last :
                domain_size_rest -= n_i
            for domain_size_cur in range(domain_size_rest+1) :
                for partition_cur in multinomial(n_cells, domain_size_cur):
                    coef = MultinomialCoefficients.coef(partition_cur)
                    val_cur = Rational(1, 1) * coef * MultinomialCoefficients.comb(domain_size_rest, domain_size_cur)
                    
                    for i in range(n_cells):
                        n_i = partition_cur[i]
                        if n_i == 0:
                            continue
                        val_cur *= cell_graph.get_cell_weight(cells[i]) ** n_i
                        val_cur *= get_two_table_weight_with_v(cell_graph, cells[i], cells[i], special_pred, v) ** (n_i * (n_i - 1) // 2)
                        for j in range(i+1, n_cells):
                            n_j = partition_cur[j]
                            if n_j == 0:
                                continue
                            val_cur *= get_two_table_weight_with_v(cell_graph, cells[i], cells[j], special_pred, v) ** (n_i * n_j)
                        for j in range(n_cells):
                            n_j = partition_last[j]
                            if n_j == 0:
                                continue
                            evidence_special_pred_false = frozenset([
                                AtomicFormula(special_pred,(Const('a'), Const('b')),False),
                                AtomicFormula(special_pred,(Const('b'), Const('a')),False) ])
                            val_cur *= cell_graph.get_two_table_weight((cells[i], cells[j]), evidence_special_pred_false) ** (n_i * n_j)
                    
                    partition_new = tuple(map(sum,zip(partition_last,partition_cur)))
                    if (dp_cur.get(partition_new)==None) :
                        dp_cur[partition_new]=Rational(0,1)
                    dp_cur[partition_new] += val_last * val_cur
        
        evaluate_u = sum( dp_cur[partition]
            for partition in multinomial(n_cells, domain_size)
                if dp_cur.get(partition)!=None )
        evaluate_u = constraint.decode_poly(evaluate_u)
        evaluate.append(evaluate_u)
        dp_last = dp_cur.copy()
        dp_cur.clear()
    
    logger.info('evaluate when v=%s: %s',v,evaluate)
    return evaluate


def F_poly(context: WFOMCContext, special_pred: Pred) -> np.poly1d :
    '''
    Return F-poly in the form of symengine expr
    '''
    domain_size = len(context.domain)
    cell_graph = CellGraph(context.formula, context.get_weight)
    MultinomialCoefficients.setup(domain_size)
    
    evaluate = F_poly_evaluate(cell_graph, domain_size, special_pred, context.cardinality_constraint)
    f_poly = lagrange.lagrange_1d(domain_size, evaluate)
    return f_poly


def Tutte_poly(context: WFOMCContext, special_pred: Pred) -> np.array :
    '''
    Return Tutte poly in the form of symengine expr
    '''
    domain_size = len(context.domain)
    cell_graph = CellGraph(context.formula, context.get_weight)
    MultinomialCoefficients.setup(domain_size)
    
    evaluate = []
    for v in range(domain_size+1) :
        evaluate.append(F_poly_evaluate(cell_graph, domain_size, special_pred, context.cardinality_constraint, v))
    evaluate = [[evaluate[v][u] for v in range(domain_size+1)] for u in range(domain_size+1)]
    tutte_poly = lagrange.lagrange_2d(domain_size, evaluate)
    return tutte_poly


def parse_args():
    parser = argparse.ArgumentParser(
        description='F-polynomial for MLN',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--output_dir', '-o', type=str,
                        default='./check-points')
    parser.add_argument('--func', '-f', type=Func,
                        choices=list(Func), default=Func.FPOLY,
                        help='the function: fpoly or tutte')
    parser.add_argument('--pred', '-p', type=str, required=True,
                        help='the special binary predicate')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if args.debug:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)
    logzero.logfile('{}/log.txt'.format(args.output_dir), mode='w')
    
    with Timer() as t:
        problem = parse_input(args.input)
    context = WFOMCContext(problem)
    logger.info('Parse input: %ss', t)
    
    if (args.func==Func.FPOLY) :
        f_poly = F_poly(context, Pred(args.pred, 2))
        logger.info('F-poly: \n%s', f_poly)
    elif (args.func==Func.TUTTE) :
        tutte_poly = Tutte_poly(context, Pred(args.pred, 2))
        logger.info('Tutte poly: \n%s', tutte_poly)