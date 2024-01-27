from symengine import Rational, var, Expr, Symbol
from .polynomial import Rational, Poly


def lagrange_1d(x, evaluate, x_symbol = Symbol('x')) :
    '''
    Input: [x1, ..., xn], [f(x1), ..., f(xn)]
        (x should be pairwise different)
    Output: f
    Default symbol for the variable: x
    '''
    n = len(x)
    res = 0
    for i in range(n) :
        L_numerator, L_denominator = 1, 1
        for j in range(n) :
            if i!=j :
                L_numerator *= x_symbol-x[j]
                L_denominator *= x[i]-x[j]
        res += L_numerator / L_denominator * evaluate[i]
    return res.expand()


def lagrange_2d(x, y, evaluate, x_symbol = Symbol('x'), y_symbol = Symbol('y')) :
    '''
    Input: [x1, ..., xn], [y1, ..., ym], [[f(x1,y1), f(x1,y2), ..., f(x1,ym)],
                                          ...
                                          [f(xn,y1), f(x1,y2), ..., f(xn,ym)]]
        (x should be pairwise different)
        (y should be pairwise different)
    Output: f
    Default symbol for the variables: x, y
    '''
    n = len(x)
    m = len(y)
    Lx_numerator, Ly_numerator, Lx_denominator, Ly_denominator = [], [], [], []
    for i in range(n) :
        numerator, denominator = 1, 1
        for j in range(n) :
            if i!=j :
                numerator *= x_symbol-x[j]
                denominator *= x[i]-x[j]
        Lx_numerator.append(numerator)
        Lx_denominator.append(denominator)
    for i in range(m) :
        numerator, denominator = 1, 1
        for j in range(m) :
            if i!=j :
                numerator *= y_symbol-y[j]
                denominator *= y[i]-y[j]
        Ly_numerator.append(numerator)
        Ly_denominator.append(denominator)
    
    res = 0
    for i in range(n) :
        for j in range(m) :
            Lx = Lx_numerator[i] / Lx_denominator[i]
            Ly = Ly_numerator[j] / Ly_denominator[j]
            res += Lx * Ly * evaluate[i][j]
    return res.expand()