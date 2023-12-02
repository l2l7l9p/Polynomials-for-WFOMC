import logging
import logzero

from logzero import logger
from symengine import Symbol, sqrt

from sampling_fo2.utils import PREDS_FOR_EXISTENTIAL, MultinomialCoefficients, multinomial, \
    multinomial_less_than, RingElement, Rational, round_rational, lagrange
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


    def _get_typed_two_table_weight(self, cell_i: Cell, cell_j: Cell, weight_type: int) :
        if (weight_type==0) :
            # ordinary two table weight
            return self.cell_graph.get_two_table_weight((cell_i, cell_j))
        elif (weight_type==-1) :
            # two table weight when no edge between cell_i and cell_j
            evidence_special_pred_false = frozenset([
                AtomicFormula(self.special_pred,(Const('a'), Const('b')),False),
                AtomicFormula(self.special_pred,(Const('b'), Const('a')),False) ])
            return self.cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_false)
        elif (weight_type==1) :
            # two table weight when cell_i <- cell_j is not allowed
            evidence_special_pred_directed = frozenset([
                AtomicFormula(self.special_pred,(Const('b'), Const('a')),False) ])
            return self.cell_graph.get_two_table_weight((cell_i, cell_j), evidence_special_pred_directed)


    def _F_poly_evaluate(self) -> list :
        cells = self.cell_graph.get_cells()
        n_cells = len(cells)
        dp_last : dict[tuple[int], Rational] = {tuple([0 for i in cells]) : 1}
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
                            if n_i == 0:
                                continue
                            val_cur *= self.cell_graph.get_cell_weight(cell_i) ** n_i
                            val_cur *= self._get_typed_two_table_weight(cell_i, cell_i, 0) ** (n_i * (n_i - 1) // 2)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_cur)):
                                if (j<=i or n_j==0):
                                    continue
                                val_cur *= self._get_typed_two_table_weight(cell_i, cell_j, 0) ** (n_i * n_j)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_last)):
                                if n_j == 0:
                                    continue
                                val_cur *= self._get_typed_two_table_weight(cell_i, cell_j, -1) ** (n_i * n_j)
                        
                        partition_new = tuple(map(sum,zip(partition_last,partition_cur)))
                        dp_cur.update({partition_new: dp_cur.get(partition_new, 0)+val_last*val_cur})
            
            evaluate_u = sum( dp_cur.get(partition, 0)
                for partition in multinomial(n_cells, self.domain_size) )
            evaluate.append(evaluate_u)
            dp_last = dp_cur.copy()
            dp_cur.clear()
        
        return evaluate


    def _G_poly_evaluate(self) :
        cells = self.cell_graph.get_cells()
        n_cells = len(cells)
        evaluate = [[] for i in range(self.domain_size+1)]
        
        dp_k_last : dict[tuple[int], Rational] = {tuple([0 for i in cells]) : 1}
        dp_k_cur : dict[tuple[int], Rational] = {}
        for k in range(self.domain_size+1) :
            # dp for k direction
            for partition_last, val_last in dp_k_last.items() :
                domain_size_last = sum(partition_last)
                domain_size_rest = self.domain_size-domain_size_last
                for domain_size_cur in range(domain_size_rest+1) :
                    for partition_cur in multinomial(n_cells, domain_size_cur):
                        coef = MultinomialCoefficients.coef(partition_cur)
                        val_cur = Rational(1, 1) * coef * MultinomialCoefficients.comb(domain_size_last+domain_size_cur, domain_size_cur)
                        
                        for i, (cell_i, n_i) in enumerate(zip(cells, partition_cur)):
                            if n_i == 0:
                                continue
                            val_cur *= self.cell_graph.get_cell_weight(cell_i) ** n_i
                            val_cur *= self._get_typed_two_table_weight(cell_i, cell_i, 0) ** (n_i * (n_i - 1) // 2)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_cur)):
                                if (j<=i or n_j==0):
                                    continue
                                val_cur *= self._get_typed_two_table_weight(cell_i, cell_j, 0) ** (n_i * n_j)
                            for j, (cell_j, n_j) in enumerate(zip(cells, partition_last)):
                                if n_j == 0:
                                    continue
                                val_cur *= self._get_typed_two_table_weight(cell_i, cell_j, 1) ** (n_i * n_j)
                        
                        partition_new = tuple(map(sum,zip(partition_last,partition_cur)))
                        dp_k_cur.update({partition_new: dp_k_cur.get(partition_new, 0)+val_last*val_cur})
            
            # dp for l direction
            dp_l_last : dict[tuple[int], Rational] = {tuple([0 for i in cells]) : 1}
            dp_l_cur : dict[tuple[int], Rational] = {}
            for l in range(self.domain_size+1) :
                for partition_last, val_last in dp_l_last.items() :
                    domain_size_last = sum(partition_last)
                    domain_size_rest = self.domain_size-domain_size_last
                    for domain_size_cur in range(domain_size_rest+1) :
                        for partition_cur in multinomial(n_cells, domain_size_cur):
                            val_cur = Rational(1, 1) * MultinomialCoefficients.comb(domain_size_last+domain_size_cur, domain_size_cur)
                            
                            val_cur *= dp_k_cur[partition_cur]
                            for i, (cell_i, n_i) in enumerate(zip(cells, partition_cur)):
                                if n_i == 0:
                                    continue
                                for j, (cell_j, n_j) in enumerate(zip(cells, partition_last)):
                                    if n_j == 0:
                                        continue
                                    val_cur *= self._get_typed_two_table_weight(cell_i, cell_j, -1) ** (n_i * n_j)
                            
                            partition_new = tuple(map(sum,zip(partition_last,partition_cur)))
                            dp_l_cur.update({partition_new: dp_l_cur.get(partition_new, 0)+val_last*val_cur})
                
                evaluate_l = sum( dp_l_cur.get(partition, 0)
                    for partition in multinomial(n_cells, self.domain_size) )
                evaluate[k].append(evaluate_l)
                dp_l_last = dp_l_cur.copy()
                dp_l_cur.clear()
            
            dp_k_last = dp_k_cur.copy()
            dp_k_cur.clear()
        
        return evaluate


    def F_poly(self) :
        '''
        Return F-poly in the form of symengine expr
        To evaluate f_poly(x0) for some x0, use f_poly.subs(Symbol('x'),x0)
        '''
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._F_poly_evaluate()
        evaluate = list(map(self.context.decode_result, evaluate))
        logger.info('evaluate: %s',evaluate)
        f_poly = lagrange.lagrange_1d(list(range(self.domain_size+1)), evaluate)
        logger.info('F-poly: %s', f_poly)
        return f_poly


    def G_poly(self) :
        '''
        Return G-poly in the form of symengine expr
        To evaluate g_poly(x0, y0) for some x0 and y0, use g_poly.subs({Symbol('x'): x0, Symbol('y'): y0})
        '''
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._G_poly_evaluate()
        evaluate = [list(map(self.context.decode_result, row)) for row in evaluate]
        logger.info('evaluate: %s',evaluate)
        g_poly = lagrange.lagrange_2d(list(range(self.domain_size+1)), list(range(self.domain_size+1)), evaluate)
        logger.info('G-poly: %s', g_poly)
        return g_poly


    def Prop2(self) :
        '''
        Return the poly in Proposition 2 in the form of symengine expr
        To evaluate prop2(x0,y0) for some x0 and y0, use prop2.subs({Symbol('x'): x0, Symbol('y'): y0})
        '''
        special_pred_original_weight = self.context.get_weight(self.special_pred)
        self.context.weights.update({self.special_pred: (special_pred_original_weight[0]*sqrt(Symbol('y')), special_pred_original_weight[1])})
        
        self.cell_graph = CellGraph(self.context.formula, self.context.get_weight)
        
        evaluate = self._F_poly_evaluate()
        evaluate = list(map(self.context.decode_result, evaluate))
        logger.info('evaluate: %s',evaluate)
        prop2 = lagrange.lagrange_1d(list(range(self.domain_size+1)), evaluate)
        logger.info('Prop2: %s', prop2)
        return prop2


    def Tutte_poly(self) :
        '''
        Return Tutte poly in the form of symengine expr
        To evaluate tutte_poly(x0,y0) for some x0 and y0, use tutte_poly.subs({Symbol('x'): x0, Symbol('y'): y0})
        '''
        prop2 = self.Prop2()
        x, y, u, v = Symbol('x'), Symbol('y'), Symbol('u'), Symbol('v')
        tutte_poly_mid = (prop2.subs({x: u*v-1, y: v}).expand() / (u*v**self.domain_size)).expand()
        tutte_poly = tutte_poly_mid.subs({u: x-1, v: y-1}).expand()
        logger.info('Tutte poly: %s', tutte_poly)
        return tutte_poly