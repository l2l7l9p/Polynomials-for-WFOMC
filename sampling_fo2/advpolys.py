import logging
import logzero

from logzero import logger
from symengine import Symbol, sqrt

from sampling_fo2.utils.multinomial import *
from sampling_fo2.utils import RingElement, Rational, round_rational, lagrange
from sampling_fo2.utils.polynomial import coeff_dict, create_vars, expand
from sampling_fo2.cell_graph import CellGraph, Cell
from sampling_fo2.context import WFOMCContext
from sampling_fo2.fol.syntax import Const, Pred, QFFormula, AtomicFormula


class AdvPolys() :
    def __init__(self, context: WFOMCContext, special_pred: Pred) :
        self.context = context
        self.special_pred = special_pred
        self.domain_size = len(context.domain)
        MultinomialCoefficients.setup(self.domain_size)


    def _dp_preprocess(self, cells, is_sscp) :
        '''
        Precompute the conditional two table weights and w_in
        '''
        self.cond_two_table_weight_ordinary = {}    # ordinary two table weight
        self.cond_two_table_weight_noedge = {}    # two table weight when no edge between cell_i and cell_j
        self.cond_two_table_weight_directed = {}    # two table weight when cell_i <- cell_j is not allowed
        evidence_special_pred_false = frozenset([
            AtomicFormula(self.special_pred,(Const('a'), Const('b')),False),
            AtomicFormula(self.special_pred,(Const('b'), Const('a')),False) ])
        evidence_special_pred_directed = frozenset([
            AtomicFormula(self.special_pred,(Const('b'), Const('a')),False) ])
        
        n_cells = len(cells)
        for cell_i in cells:
            for cell_j in cells:
                self.cond_two_table_weight_ordinary[(cell_i, cell_j)] = self.cell_graph.get_two_table_weight((cell_i, cell_j))
                self.cond_two_table_weight_noedge[(cell_i, cell_j)] = self.cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_false)
                self.cond_two_table_weight_directed[(cell_i, cell_j)] = self.cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_directed)
        
        self.w_in = {}
        cond_two_table_weight = self.cond_two_table_weight_noedge if is_sscp else self.cond_two_table_weight_ordinary
        for partition in multinomial_less_than(n_cells, self.domain_size) :
            val = Rational(1, 1) * MultinomialCoefficients.coef(partition)
            for i, (cell_i, n_i) in enumerate(zip(cells, partition)):
                val *= self.cell_graph.get_cell_weight(cell_i) ** n_i
                val *= cond_two_table_weight[(cell_i, cell_i)] ** (n_i * (n_i - 1) // 2)
                for j, (cell_j, n_j) in enumerate(zip(cells, partition)):
                    if (j>=i) : break
                    val *= cond_two_table_weight[(cell_i, cell_j)] ** (n_i * n_j)
            self.w_in[partition] = val


    def _WCP_evaluate(self) -> list :
        cells = self.cell_graph.get_cells()
        n_cells = len(cells)
        self._dp_preprocess(cells, False)
        evaluate = []
        for u in range(self.domain_size+1) :
            if u==0 :
                dp = {partition: self.w_in[partition] for partition in multinomial_less_than(n_cells, self.domain_size)}
            else :
                for size_new in range(self.domain_size,-1,-1) :
                    for partition_new in multinomial(n_cells, size_new) :
                        val = 0
                        for partition_cur in multinomial_less_than_multinomial(n_cells, partition_new) :
                            size_cur = sum(partition_cur)
                            partition_last = tuple((partition_new[i]-partition_cur[i] for i in range(n_cells)))
                            val_cur = MultinomialCoefficients.comb(size_new, size_cur) * self.w_in[partition_cur]
                            for cell_i, n_i in zip(cells, partition_cur):
                                for cell_j, n_j in zip(cells, partition_last):
                                    val_cur *= self.cond_two_table_weight_noedge[(cell_i, cell_j)] ** (n_i * n_j)
                            val += dp[partition_last]*val_cur
                        dp[partition_new] = val
            
            evaluate_u = sum(dp[partition] for partition in multinomial(n_cells, self.domain_size))
            evaluate.append(evaluate_u)
        
        self.w_in = None
        return evaluate


    def _SCP_evaluate(self, is_sscp) :
        cells = self.cell_graph.get_cells() if not is_sscp else self.cell_graph.get_cells(lambda x: not x.is_positive(self.special_pred))
        n_cells = len(cells)
        self._dp_preprocess(cells, is_sscp)
        evaluate = [[] for i in range(self.domain_size+1)]
        
        for v in range(self.domain_size+1) :
            # dp for v direction
            if v==0 :
                dp_v = {partition: self.w_in[partition] for partition in multinomial_less_than(n_cells, self.domain_size)}
            else :
                for size_new in range(self.domain_size,-1,-1) :
                    for partition_new in multinomial(n_cells, size_new) :
                        val = 0
                        for partition_cur in multinomial_less_than_multinomial(n_cells, partition_new) :
                            size_cur = sum(partition_cur)
                            partition_last = tuple((partition_new[i]-partition_cur[i] for i in range(n_cells)))
                            val_cur = MultinomialCoefficients.comb(size_new, size_cur) * self.w_in[partition_cur]
                            for cell_i, n_i in zip(cells, partition_cur):
                                for cell_j, n_j in zip(cells, partition_last):
                                    val_cur *= self.cond_two_table_weight_directed[(cell_i, cell_j)] ** (n_i * n_j)
                            val += dp_v[partition_last]*val_cur
                        dp_v[partition_new] = val
            
            # dp for u direction
            for u in range(self.domain_size+1) :
                if u==0 :
                    dp_u = {partition: dp_v[partition] for partition in multinomial_less_than(n_cells, self.domain_size)}
                else :
                    for size_new in range(self.domain_size,-1,-1) :
                        for partition_new in multinomial(n_cells, size_new) :
                            val = 0
                            for partition_cur in multinomial_less_than_multinomial(n_cells, partition_new) :
                                size_cur = sum(partition_cur)
                                partition_last = tuple((partition_new[i]-partition_cur[i] for i in range(n_cells)))
                                val_cur = MultinomialCoefficients.comb(size_new, size_cur) * dp_v[partition_cur]
                                for cell_i, n_i in zip(cells, partition_cur):
                                    for cell_j, n_j in zip(cells, partition_last):
                                        val_cur *= self.cond_two_table_weight_noedge[(cell_i, cell_j)] ** (n_i * n_j)
                                val += dp_u[partition_last]*val_cur
                            dp_u[partition_new] = val
                
                evaluate_u = sum(dp_u[partition] for partition in multinomial(n_cells, self.domain_size))
                evaluate[v].append(evaluate_u)
        
        evaluate = list(map(list, zip(*evaluate)))
        self.w_in = None
        return evaluate


    def WCP(self) :
        '''
        Return the Weak Connectedness Polynomial in the form of symengine expr
        To evaluate wcp(u0) for some u0, use wcp.subs(Symbol('u'),u0)
        '''
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._WCP_evaluate()
        evaluate = list(map(self.context.decode_result, evaluate))
        logger.info('evaluate: %s',evaluate)
        wcp = lagrange.lagrange_1d(list(range(self.domain_size+1)), evaluate, Symbol('u'))
        logger.info('WCP: %s', wcp)
        return wcp


    def SCP(self, strictness) :
        '''
        Return the Strong Connectedness Polynomial in the form of symengine expr
        To evaluate scp(u0, v0) for some u0 and v0, use scp.subs({Symbol('u'): u0, Symbol('v'): v0})
        '''
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._SCP_evaluate(strictness)
        evaluate = [list(map(self.context.decode_result, row)) for row in evaluate]
        logger.info('evaluate: %s',evaluate)
        scp = lagrange.lagrange_2d(list(range(self.domain_size+1)), list(range(self.domain_size+1)), evaluate, Symbol('u'), Symbol('v'))
        logger.info('SCP: %s', scp)
        return scp


    def extended_WCP(self) :
        '''
        Return the extended Weak Connectedness Polynomial in the form of symengine expr
        To evaluate ewcp(u0, v0) for some u0 and v0, use ewcp.subs({Symbol('u'): u0, Symbol('v'): v0})
        '''
        special_pred_original_weight = self.context.get_weight(self.special_pred)
        self.context.weights.update({self.special_pred: (special_pred_original_weight[0]*Symbol('v'), special_pred_original_weight[1])})
        
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._WCP_evaluate()
        evaluate = list(map(self.context.decode_result, evaluate))
        logger.info('evaluate: %s',evaluate)
        ewcp = lagrange.lagrange_1d(list(range(self.domain_size+1)), evaluate, Symbol('u'))
        logger.info('extended WCP: %s', ewcp)
        self.context.weights.update({self.special_pred: special_pred_original_weight})
        return ewcp


    def Tutte_poly(self) :
        '''
        Return the Tutte Polynomial in the form of symengine expr
        To evaluate tutte_poly(x0,y0) for some x0 and y0, use tutte_poly.subs({Symbol('x'): x0, Symbol('y'): y0})
        '''
        ewcp = self.extended_WCP()
        x, y, u, v = Symbol('x'), Symbol('y'), Symbol('u'), Symbol('v')
        tutte_poly_mid = (ewcp.subs({u: x*y-1, v: sqrt(y)}).expand() / (x*y**self.domain_size)).expand()
        tutte_poly = tutte_poly_mid.subs({x: x-1, y: y-1}).expand()
        logger.info('Tutte poly: %s', tutte_poly)
        return tutte_poly