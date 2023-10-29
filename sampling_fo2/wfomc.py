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

from sampling_fo2.cell_graph.cell_graph import OptimizedCellGraph

from sampling_fo2.utils import PREDS_FOR_EXISTENTIAL, MultinomialCoefficients, multinomial, \
    multinomial_less_than, RingElement, Rational, round_rational
from sampling_fo2.cell_graph import CellGraph, Cell
from sampling_fo2.context import WFOMCContext
from sampling_fo2.parser import parse_input
from sampling_fo2.fol.syntax import Const, Pred, QFFormula
from sampling_fo2.utils.polynomial import coeff_dict, create_vars, expand


class Algo(Enum):
    STANDARD = 'standard'
    FASTER = 'faster'

    def __str__(self):
        return self.value


def get_config_weight_standard_faster(config: list[int],
                                      cell_weights: list[RingElement],
                                      edge_weights: list[list[RingElement]]) \
        -> RingElement:
    res = Rational(1, 1)
    for i, n_i in enumerate(config):
        if n_i == 0:
            continue
        n_i = Rational(n_i, 1)
        res *= cell_weights[i] ** n_i
        res *= edge_weights[i][i] ** (n_i * (n_i - 1) // Rational(2, 1))
        for j, n_j in enumerate(config):
            if j <= i:
                continue
            if n_j == 0:
                continue
            n_j = Rational(n_j, 1)
            res *= edge_weights[i][j] ** (n_i * n_j)
    return res


def get_config_weight_standard(cell_graph: CellGraph,
                               cell_config: dict[Cell, int]) -> RingElement:
    res = Rational(1, 1)
    for cell, n in cell_config.items():
        if n > 0:
            # NOTE: nullary weight is multiplied once
            res = res * cell_graph.get_nullary_weight(cell)
            break
    for i, (cell_i, n_i) in enumerate(cell_config.items()):
        if n_i == 0:
            continue
        res = res * cell_graph.get_cell_weight(cell_i) ** n_i
        res = res * cell_graph.get_two_table_weight(
            (cell_i, cell_i)
        ) ** (n_i * (n_i - 1) // 2)
        for j, (cell_j, n_j) in enumerate(cell_config.items()):
            if j <= i:
                continue
            if n_j == 0:
                continue
            res = res * cell_graph.get_two_table_weight(
                (cell_i, cell_j)
            ) ** (n_i * n_j)
    # logger.debug('Config weight: %s', res)
    return res


def faster_wfomc(formula: QFFormula,
               domain: set[Const],
               get_weight: Callable[[Pred], tuple[RingElement, RingElement]]) -> RingElement:
    domain_size = len(domain)
    with Timer() as t:
        opt_cell_graph = OptimizedCellGraph(formula, get_weight, domain_size)
    logger.info('Optimized cell graph construction time: %s', t.elapsed)

    cliques = opt_cell_graph.cliques
    nonind = opt_cell_graph.nonind
    i2_ind = opt_cell_graph.i2_ind
    nonind_map = opt_cell_graph.nonind_map

    res = Rational(0, 1)
    with Timer() as t:
        for partition in multinomial_less_than(len(nonind), domain_size):
            mu = tuple(partition)
            if sum(partition) < domain_size:
                mu = mu + (domain_size - sum(partition),)
            coef = MultinomialCoefficients.coef(mu)
            body = Rational(1, 1)

            for i, clique1 in enumerate(cliques):
                for j, clique2 in enumerate(cliques):
                    if i in nonind and j in nonind:
                        if i < j:
                            body = body * opt_cell_graph.get_two_table_weight(
                                (clique1[0], clique2[0])
                            ) ** (partition[nonind_map[i]] *
                                  partition[nonind_map[j]])

            for l in nonind:
                body = body * opt_cell_graph.get_J_term(
                    l, partition[nonind_map[l]]
                )
                body = body * opt_cell_graph.get_cell_weight(
                    cliques[l][0]
                ) ** partition[nonind_map[l]]

            opt_cell_graph.setup_term_cache()
            mul = opt_cell_graph.get_term(len(i2_ind), 0, partition)
            res = res + coef * mul * body
    logger.info('WFOMC time: %s', t.elapsed)

    return res


def standard_wfomc(formula: QFFormula,
                   domain: set[Const],
                   get_weight: Callable[[Pred], tuple[RingElement, RingElement]]) -> RingElement:
    cell_graph = CellGraph(formula, get_weight)
    # cell_graph.show()
    cells = cell_graph.get_cells()
    n_cells = len(cells)
    domain_size = len(domain)
    MultinomialCoefficients.setup(domain_size)

    res = Rational(0, 1)
    for partition in multinomial(n_cells, domain_size):
        coef = MultinomialCoefficients.coef(partition)
        cell_config = dict(zip(cells, partition))
        # logger.debug(
        #     '=' * 15 + ' Compute WFOMC for the partition %s ' + '=' * 15,
        #     dict(filter(lambda x: x[1] != 0, cell_config.items())
        # ))
        res = res + coef * get_config_weight_standard(
            cell_graph, cell_config
        )
    return res


def precompute_ext_weight(cell_graph: CellGraph, domain_size: int,
                          context: WFOMCContext) \
        -> dict[frozenset[tuple[Cell, frozenset[Pred], int]], RingElement]:
    existential_weights = defaultdict(lambda: Rational(0, 1))
    cells = cell_graph.get_cells()
    eu_configs = []
    for cell in cells:
        config = []
        for domain_pred, tseitin_preds in context.domain_to_evidence_preds.items():
            if cell.is_positive(domain_pred):
                config.append((
                    cell.drop_preds(
                        prefixes=PREDS_FOR_EXISTENTIAL), tseitin_preds
                ))
        eu_configs.append(config)

    cell_weights, edge_weights = cell_graph.get_all_weights()

    for partition in multinomial(len(cells), domain_size):
        # res = get_config_weight_standard(
        #     cell_graph, dict(zip(cells, partition))
        # )
        res = get_config_weight_standard_faster(
            partition, cell_weights, edge_weights
        )
        eu_config = defaultdict(lambda: 0)
        for idx, n in enumerate(partition):
            for config in eu_configs[idx]:
                eu_config[config] += n
        eu_config = dict(
            (k, v) for k, v in eu_config.items() if v > 0
        )
        existential_weights[
            frozenset((*k, v) for k, v in eu_config.items())
        ] += (Rational(MultinomialCoefficients.coef(partition), 1) * res)
    # remove duplications
    for eu_config in existential_weights.keys():
        dup_factor = Rational(MultinomialCoefficients.coef(
            tuple(c[2] for c in eu_config)
        ), 1)
        existential_weights[eu_config] /= dup_factor
    return existential_weights


def wfomc(context: WFOMCContext, algo: Algo = Algo.STANDARD) -> Rational:
    if algo == Algo.STANDARD:
        res = standard_wfomc(
            context.formula, context.domain, context.get_weight
        )
    elif algo == Algo.FASTER:
        res = faster_wfomc(
            context.formula, context.domain, context.get_weight
        )
    return context.decode_result(res)


def count_distribution(context: WFOMCContext, preds: list[Pred],
                       algo: Algo = Algo.STANDARD) \
        -> dict[tuple[int, ...], Rational]:
    pred2weight = {}
    pred2sym = {}
    syms = create_vars('x0:{}'.format(len(preds)))
    for sym, pred in zip(syms, preds):
        if pred in pred2weight:
            continue
        weight = context.get_weight(pred)
        pred2weight[pred] = (weight[0] * sym, weight[1])
        pred2sym[pred] = sym
    context.weights.update(pred2weight)
    if algo == Algo.STANDARD:
        res = standard_wfomc(
            context.formula, context.domain, context.get_weight
        )
    elif algo == Algo.FASTER:
        res = faster_wfomc(
            context.formula, context.domain, context.get_weight
        )

    symbols = [pred2sym[pred] for pred in preds]
    count_dist = {}
    res = expand(res)
    for degrees, coef in coeff_dict(res, symbols):
        count_dist[degrees] = coef
    return count_dist


def parse_args():
    parser = argparse.ArgumentParser(
        description='WFOMC for MLN',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--output_dir', '-o', type=str,
                        default='./check-points')
    parser.add_argument('--algo', '-a', type=Algo,
                        choices=list(Algo), default=Algo.FASTER)
    parser.add_argument('--domain_recursive',
                        action='store_true', default=False,
                        help='use domain recursive algorithm '
                             '(only for existential quantified MLN)')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # import sys
    # sys.setrecursionlimit(int(1e6))
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

    res = wfomc(
        context, algo=args.algo
    )
    logger.info('WFOMC (arbitrary precision): %s', res)
    round_val = round_rational(res)
    logger.info('WFOMC (round): %s (exp(%s))', round_val, round_val.ln())
