from lark import Lark
from sampling_fo2.network.mln import MLN
from sampling_fo2.fol.syntax import Const, Pred
from sampling_fo2.network.constraint import CardinalityConstraint
from sampling_fo2.parser.mln_grammar import grammar
from sampling_fo2.utils import Rational


from sampling_fo2.parser.fol_parser import FOLTransformer
from sampling_fo2.problems import MLNProblem

from sampling_fo2.fol.syntax import *

class MLNTransformer(FOLTransformer):

    def domain_elements(self, args):
        return list(args)

    def int_domain(self, args):
        return int(args[0])

    def element(self, args):
        return args[0].value

    def set_domain(self, args):
        return set(args[0])

    def domain_name(self, args):
        return args[0].value

    def domain(self, args):
        domain_name, domain_spec = args
        if isinstance(domain_spec, int):
            domain_spec = set(f'{domain_name}{i}' for i in range(domain_spec))
        return (domain_name, domain_spec)

    def cardinality_param(self, args):
        return int(args[0])

    def cardinality(self, args):
        pred, op, param = args
        return (pred, op, param)

    def cardinalities(self, args):
        if len(args) == 0:
            return None
        return list(args)

    def weighting(self, args):
        return float(args[0])

    def rules(self, args):
        rules = args
        weightings = []
        formulas = []
        for w, r in rules:
            weightings.append(w)
            formulas.append(r)
        return weightings, formulas

    def rule(self, args):
        w, r = args[0]
        return w, r

    def hard_rule(self, args):
        return float('inf'), args[0]

    def soft_rule(self, args):
        return args[0], args[1]

    def mln(self, args):
        rules = args[0]
        domain = args[1][1] # Only one definition domain is supported
        cardinalities = args[2]

        cc_constraints: dict[Pred, tuple[str, int]] = dict()
        if cardinalities is not None:
            for pred, op, param in cardinalities:
                pred = self.name2pred.get(pred, None)
                if not pred:
                    raise ValueError(f'Predicate {pred} not found')
                cc_constraints[pred] = (op, param)
            cardinality_constraint = CardinalityConstraint(cc_constraints)
        else:
            cardinality_constraint = None

        return rules, domain, cardinality_constraint

def parse(text: str) -> MLNProblem:
    mln_parser = Lark(grammar,
                        start='mln')
    tree = mln_parser.parse(text)
    (rules, domain, cardinality_constraint) = MLNTransformer().transform(tree)

    return MLNProblem(
        rules,
        domain,
        cardinality_constraint
    )
