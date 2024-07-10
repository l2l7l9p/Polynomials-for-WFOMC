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
    WCP = 'wcp'
    NSCP = 'nscp'
    SSCP = 'sscp'
    EWCP = 'ewcp'
    TUTTE = 'tutte'

    def __str__(self):
        return self.value


def parse_args():
    parser = argparse.ArgumentParser(
        description='Computing Polynomials from a C2 sentence',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--output_dir', '-o', type=str,
                        default='./check-points')
    parser.add_argument('--func', '-f', type=Func,
                        choices=list(Func), default=Func.WCP,
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
    
    problem = parse_input(args.input)
    context = WFOMCContext(problem)
    polys = AdvPolys(context, Pred(args.pred, 2))
    
    if (args.func==Func.WCP) :
        with Timer() as t:
            wcp = polys.WCP()
        logger.info('Time of computing WCP: %ss', t)
        
    elif (args.func==Func.NSCP) :
        with Timer() as t:
            scp = polys.SCP(False)
        logger.info('Time of computing NSCP: %ss', t)
        # scp = scp.subs({Symbol('u'): 0, Symbol('v'): -2}).expand()
        # logger.info(scp)
    
    elif (args.func==Func.SSCP) :
        with Timer() as t:
            scp = polys.SCP(True)
        logger.info('Time of computing SSCP: %ss', t)
        scp = scp.subs({Symbol('u'): 0, Symbol('v'): -2}).expand()
        logger.info(scp)
    
    elif (args.func==Func.EWCP) :
        with Timer() as t:
            ewcp = polys.extended_WCP()
        logger.info('Time of computing extended WCP: %ss', t)
    
    elif (args.func==Func.TUTTE) :
        tutte_poly = polys.Tutte_poly()
        logger.info('Tutte poly at (1,1): %s', tutte_poly.subs({Symbol('x'): 1, Symbol('y'): 1}))