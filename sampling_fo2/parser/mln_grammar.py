from .fol_grammar import function_free_logic_grammar

domain_grammar = r"""
    domain: domain_name "=" domain_spec
    domain_name: CNAME
    ?domain_spec: INT               -> int_domain
        | ("{" domain_elements "}") -> set_domain
    domain_elements: element ("," element)*
    element: CNAME
"""

cardinality_grammar = r"""
    cardinalities: cardinality*
    cardinality: "|" predicate "|" comparator cardinality_param
    cardinality_param: INT
"""

rule_grammar = r"""
    rules: rule*
    rule: hard_rule | soft_rule
    hard_rule: ffl "."
    soft_rule: weighting ffl
    weighting: SIGNED_NUMBER
""" + function_free_logic_grammar

grammar = r"""
    ?mln: rules domain cardinalities
""" + rule_grammar + domain_grammar + cardinality_grammar
