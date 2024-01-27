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


grammar = r"""
    ?wfomcs: ffl domain weightings cardinalities
    weightings: weighting*
    weighting: weight weight predicate

    weight: FLOAT | INT
""" + domain_grammar + cardinality_grammar + function_free_logic_grammar
