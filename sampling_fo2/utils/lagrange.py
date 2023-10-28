from symengine import Rational, var, Expr, Symbol
from .polynomial import Rational, Poly


def lagrange_1d(n, evaluate) :
    '''
    Input: [f(0), f(1), ..., f(n)]
    Output: f
    '''
    res = 0
    for i in range(n+1) :
        L_numerator = 1
        L_denominator = 1
        for j in range(n+1) :
            if i!=j :
                L_numerator *= Symbol('x')-j
                L_denominator *= i-j
        res += L_numerator / L_denominator * evaluate[i]
    return res.expand()


def lagrange_2d(n, evaluate) :
    '''
    Input: [[f(0,0), f(0,1), ..., f(0,n)],
            ...
            [f(n,0), f(n,1), ..., f(n,n)]]
    Output: f
    '''
    Lx_numerator = []
    Ly_numerator = []
    L_denominator = []
    for i in range(n+1) :
        numerator_x = 1
        numerator_y = 1
        denominator = 1
        for j in range(n+1) :
            if i!=j :
                numerator_x *= Symbol('x')-j
                numerator_y *= Symbol('y')-j
                denominator *= i-j
        Lx_numerator.append(numerator_x)
        Ly_numerator.append(numerator_y)
        L_denominator.append(denominator)
    
    res = 0
    for i in range(n+1) :
        for j in range(n+1) :
            Lx = Lx_numerator[i] / L_denominator[i]
            Ly = Ly_numerator[j] / L_denominator[j]
            res += Lx * Ly * evaluate[i][j]
    return res.expand()