import logging
import logzero

from logzero import logger
from symengine import Symbol, sqrt

from sampling_fo2.utils import MultinomialCoefficients, multinomial, RingElement, Rational, round_rational, lagrange
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


    def _get_conditional_two_table(self) :
        '''
        Precompute the conditional two table weights
        '''
        self.cond_two_table_weight_ordinary = {}    # ordinary two table weight
        self.cond_two_table_weight_noedge = {}    # two table weight when no edge between cell_i and cell_j
        self.cond_two_table_weight_directed = {}    # two table weight when cell_i <- cell_j is not allowed
        evidence_special_pred_false = frozenset([
            AtomicFormula(self.special_pred,(Const('a'), Const('b')),False),
            AtomicFormula(self.special_pred,(Const('b'), Const('a')),False) ])
        evidence_special_pred_directed = frozenset([
            AtomicFormula(self.special_pred,(Const('b'), Const('a')),False) ])
        
        cells = self.cell_graph.get_cells()
        for cell_i in cells:
            for cell_j in cells:
                self.cond_two_table_weight_ordinary[(cell_i, cell_j)] = self.cell_graph.get_two_table_weight((cell_i, cell_j))
                self.cond_two_table_weight_noedge[(cell_i, cell_j)] = self.cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_false)
                self.cond_two_table_weight_directed[(cell_i, cell_j)] = self.cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_directed)


    def _WCP_evaluate(self) -> list :
        cells = self.cell_graph.get_cells()
        n_cells = len(cells)
        self._get_conditional_two_table()
        dp_last : dict[tuple[int], Rational] = {tuple([0]*n_cells) : 1}
        dp_cur : dict[tuple[int], Rational] = {}
        evaluate = []
        for u in range(self.domain_size+1) :
            for partition_last, val_last in dp_last.items() :
                domain_size_last = sum(partition_last)
                domain_size_rest = self.domain_size-domain_size_last
                for domain_size_cur in range(domain_size_rest+1) :
                    for partition_cur in multinomial(n_cells, domain_size_cur):
                        coef = MultinomialCoefficients.coef(partition_cur)
                        val_cur = Rational(1, 1) * coef * MultinomialCoefficients.comb(domain_size_last+domain_size_cur, domain_size_cur)
                        
                        for i, (cell_i, n_i) in enumerate(zip(cells, partition_cur)):
                            val_cur *= self.cell_graph.get_cell_weight(cell_i) ** n_i
                            val_cur *= self.cond_two_table_weight_ordinary[(cell_i, cell_i)] ** (n_i * (n_i - 1) // 2)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_cur)):
                                if (j>i) :
                                    val_cur *= self.cond_two_table_weight_ordinary[(cell_i, cell_j)] ** (n_i * n_j)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_last)):
                                val_cur *= self.cond_two_table_weight_noedge[(cell_i, cell_j)] ** (n_i * n_j)
                        
                        partition_new = tuple(map(sum,zip(partition_last,partition_cur)))
                        dp_cur.update({partition_new: dp_cur.get(partition_new, 0)+val_last*val_cur})
            
            evaluate_u = sum( dp_cur.get(partition, 0)
                for partition in multinomial(n_cells, self.domain_size) )
            evaluate.append(evaluate_u)
            dp_last, dp_cur = dp_cur, dp_last
            dp_cur.clear()
        
        return evaluate


    def _SCP_evaluate(self) :
        cells = self.cell_graph.get_cells()
        n_cells = len(cells)
        self._get_conditional_two_table()
        evaluate = [[] for i in range(self.domain_size+1)]
        
        dp_v_last : dict[tuple[int], Rational] = {tuple([0]*n_cells) : 1}
        dp_v_cur : dict[tuple[int], Rational] = {}
        for v in range(self.domain_size+1) :
            # dp for v direction
            for partition_last, val_last in dp_v_last.items() :
                domain_size_last = sum(partition_last)
                domain_size_rest = self.domain_size-domain_size_last
                for domain_size_cur in range(domain_size_rest+1) :
                    for partition_cur in multinomial(n_cells, domain_size_cur):
                        coef = MultinomialCoefficients.coef(partition_cur)
                        val_cur = Rational(1, 1) * coef * MultinomialCoefficients.comb(domain_size_last+domain_size_cur, domain_size_cur)
                        
                        for i, (cell_i, n_i) in enumerate(zip(cells, partition_cur)):
                            val_cur *= self.cell_graph.get_cell_weight(cell_i) ** n_i
                            val_cur *= self.cond_two_table_weight_ordinary[(cell_i, cell_i)] ** (n_i * (n_i - 1) // 2)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_cur)):
                                if (j>i) :
                                    val_cur *= self.cond_two_table_weight_ordinary[(cell_i, cell_j)] ** (n_i * n_j)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_last)):
                                val_cur *= self.cond_two_table_weight_directed[(cell_i, cell_j)] ** (n_i * n_j)
                        
                        partition_new = tuple(map(sum,zip(partition_last,partition_cur)))
                        dp_v_cur.update({partition_new: dp_v_cur.get(partition_new, 0)+val_last*val_cur})
            
            # dp for u direction
            dp_u_last : dict[tuple[int], Rational] = {tuple([0]*n_cells) : 1}
            dp_u_cur : dict[tuple[int], Rational] = {}
            for u in range(self.domain_size+1) :
                for partition_last, val_last in dp_u_last.items() :
                    domain_size_last = sum(partition_last)
                    domain_size_rest = self.domain_size-domain_size_last
                    for domain_size_cur in range(domain_size_rest+1) :
                        for partition_cur in multinomial(n_cells, domain_size_cur):
                            val_cur = Rational(1, 1) * MultinomialCoefficients.comb(domain_size_last+domain_size_cur, domain_size_cur)
                            
                            val_cur *= dp_v_cur[partition_cur]
                            for i, (cell_i, n_i) in enumerate(zip(cells, partition_cur)):
                                for j, (cell_j, n_j) in enumerate(zip(cells, partition_last)):
                                    val_cur *= self.cond_two_table_weight_noedge[(cell_i, cell_j)] ** (n_i * n_j)
                            
                            partition_new = tuple(map(sum,zip(partition_last,partition_cur)))
                            dp_u_cur.update({partition_new: dp_u_cur.get(partition_new, 0)+val_last*val_cur})
                
                evaluate_u = sum( dp_u_cur.get(partition, 0)
                    for partition in multinomial(n_cells, self.domain_size) )
                evaluate[v].append(evaluate_u)
                dp_u_last, dp_u_cur = dp_u_cur, dp_u_last
                dp_u_cur.clear()
            
            dp_v_last, dp_v_cur = dp_v_cur, dp_v_last
            dp_v_cur.clear()
            # logger.info('evaluate at v=%s, %s', v, evaluate[v])
        
        evaluate = list(map(list, zip(*evaluate)))
        return evaluate


    def WCP(self) :
        '''
        Return the Weak Connectedness Polynomial in the form of symengine expr
        To evaluate wcp(x0) for some x0, use wcp.subs(Symbol('x'),x0)
        '''
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._WCP_evaluate()
        evaluate = list(map(self.context.decode_result, evaluate))
        logger.info('evaluate: %s',evaluate)
        wcp = lagrange.lagrange_1d(list(range(self.domain_size+1)), evaluate, Symbol('u'))
        logger.info('WCP: %s', wcp)
        return wcp


    def SCP(self) :
        '''
        Return the Strong Connectedness Polynomial in the form of symengine expr
        To evaluate scp(x0, y0) for some x0 and y0, use scp.subs({Symbol('x'): x0, Symbol('y'): y0})
        '''
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._SCP_evaluate()
        evaluate = [list(map(self.context.decode_result, row)) for row in evaluate]
        logger.info('evaluate: %s',evaluate)
        scp = lagrange.lagrange_2d(list(range(self.domain_size+1)), list(range(self.domain_size+1)), evaluate, Symbol('u'), Symbol('v'))
        logger.info('SCP: %s', scp)
        return scp


    def extended_WCP(self) :
        '''
        Return the extended Weak Connectedness Polynomial in the form of symengine expr
        To evaluate ewcp(x0, y0) for some x0 and y0, use ewcp.subs({Symbol('x'): x0, Symbol('y'): y0})
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