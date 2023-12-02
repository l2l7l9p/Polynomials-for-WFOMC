'''
Input: an FO formula
Output: F-poly, G-poly, Prop.2 poly or Tutte poly of it
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
from symengine import Symbol

from sampling_fo2.context import WFOMCContext
from sampling_fo2.parser import parse_input
from sampling_fo2.fol.syntax import Const, Pred, QFFormula, AtomicFormula
from sampling_fo2.utils.polynomial import coeff_dict, create_vars, expand
from sampling_fo2.advpolys import AdvPolys


class Func(Enum):
    FPOLY = 'fpoly'
    GPOLY = 'gpoly'
    TUTTE = 'tutte'
    PROP2 = 'prop2'

    def __str__(self):
        return self.value


def parse_args():
    parser = argparse.ArgumentParser(
        description='Polynomials for MLN',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--output_dir', '-o', type=str,
                        default='./check-points')
    parser.add_argument('--func', '-f', type=Func,
                        choices=list(Func), default=Func.FPOLY,
                        help='the function wanted')
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
    
    polys = AdvPolys(context, Pred(args.pred, 2))
    
    if (args.func==Func.FPOLY) :
        f_poly = polys.F_poly()
        
    elif (args.func==Func.GPOLY) :
        g_poly = polys.G_poly()
        prop3 = g_poly.subs({Symbol('x'): -2, Symbol('y'): Symbol('y')-1}).expand()
        # logger.info('Prop3: %s', prop3)
        logger.info('Num of strongly connected: %s', -prop3.as_coefficients_dict()[Symbol('y')])
        
    elif (args.func==Func.PROP2) :
        prop2 = polys.Prop2()
        
    elif (args.func==Func.TUTTE) :
        tutte_poly = polys.Tutte_poly()
        logger.info('Tutte poly at (1,1): %s', tutte_poly.subs({Symbol('x'): 1, Symbol('y'): 1}))