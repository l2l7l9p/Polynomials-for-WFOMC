from __future__ import annotations

from fractions import Fraction

from lark import Lark
from sampling_fo2.fol.sc2 import SC2, to_sc2
from sampling_fo2.fol.syntax import Const, Pred
from sampling_fo2.network.constraint import CardinalityConstraint

from sampling_fo2.parser.fol_parser import FOLTransformer
from sampling_fo2.parser.wfomcs_grammar import grammar
from sampling_fo2.problems import WFOMCSProblem
from sampling_fo2.utils import Rational


class WFOMSTransformer(FOLTransformer):

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

    def weightings(self, args):
        return dict(args)

    def weighting(self, args):
        return (args[2], (
            Rational(
                Fraction(args[0]).numerator,
                Fraction(args[0]).denominator
            ),
            Rational(
                Fraction(args[1]).numerator,
                Fraction(args[1]).denominator
            )
        ))

    def weight(self, args):
        return str(args[0])

    def wfomcs(self, args):
        sentence = args[0]
        domain = args[1][1]
        weightings = args[2]
        cardinalities = args[3]
        sentence = to_sc2(sentence)

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

        return sentence, domain, weightings, cardinality_constraint


def parse(text: str) -> \
        tuple[SC2, set[Const], dict[Pred, tuple[Rational, Rational]], CardinalityConstraint]:
    """
    Parse wfoms text into WFOMSContext
    """
    wfomcs_parser = Lark(grammar,
                        start='wfomcs')
    tree = wfomcs_parser.parse(text)
    (
        sentence,
        domain,
        weightings,
        cardinality_constraint
    ) = WFOMSTransformer().transform(tree)
    pred_weightings = dict(
        (sentence.pred_by_name(pred), weights)
        for pred, weights in weightings.items()
    )
    return WFOMCSProblem(
        sentence,
        domain,
        pred_weightings,
        cardinality_constraint
    )


if __name__ == '__main__':
    wfoms = parse(r'''
\forall X: (\forall Y: (E(X,Y) -> E(Y,X))) &
\forall X: (~E(X,X)) &
\forall X: (\exists Y: (E(X,Y)))

vertices = 10
0.1 1 E
0.2 2 F
0.3 3 G
    ''')
