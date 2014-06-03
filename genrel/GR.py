"""
genrel package for GR calculations
David Clark, Kai Smith, David Cyncynates
Case Western Reserve university
2014
"""

import numpy as np
import sympy

def christoffel(metric, metric_key):
    #symols will be a rank 3 tensor. The first index will correspond to the upper
    #index and the next two will correspond to the lower indecies.
    symbols = np.empty((4, 4, 4), dtype = type(sympy.Symbol('Test')*1))
    inverse = sympy.Matrix(metric).inv()
    for alpha in range(4):
        for beta in range(4):
            for gamma in range(4):
                total = 0
                for delta in range(4):                    
                    total += inverse[alpha, delta] * (sympy.diff(metric[delta][beta], metric_key[gamma]) + 
                            sympy.diff(metric[delta][gamma], metric_key[beta]) - sympy.diff(metric[beta][gamma], metric_key[delta]))
                symbols[alpha][beta][gamma] = sympy.cancel(total/2)
    return symbols
                    

def reimann_tensor(chris_sym, metric_key):
    #reimann tensor calculated with the first index as an upper index and all of the other ones as lower indecies
    reimann = np.empty((4, 4, 4, 4), dtype = type(sympy.Symbol('Test')*1))
    for alpha in range(4):
        for beta in range(4):
            for gamma in range(4):
                for delta in range(4):
                    total = 0
                    total += sympy.diff(chris_sym[alpha][beta][delta], metric_key[gamma])
                    total -= sympy.diff(chris_sym[alpha][beta][gamma], metric_key[delta])
                    for epsilon in range(4):
                        total += chris_sym[alpha][gamma][epsilon]*chris_sym[epsilon][beta][delta]
                        total -= chris_sym[alpha][delta][epsilon]*chris_sym[epsilon][beta][gamma]
                    reimann[alpha][beta][gamma][delta] = sympy.cancel(total)
    return reimann


def ricci_tensor(reimann):
    ricci = np.empty((4, 4), dtype = type(sympy.Symbol('Test')*1))
    for alpha in range(4):
        for beta in range(4):
            total = 0
            for gamma in range(4):
                total += reimann[gamma][alpha][gamma][beta]
            ricci[alpha][beta] = sympy.cancel(total)
    return ricci

def ricci_scalar(ricci_t, metric):
    scalar = 0
    inverse = sympy.Matrix(metric).inv()
    for alpha in range(4):
        for beta in range(4):
            scalar += inverse[alpha, beta] * ricci_t[alpha][beta]
    return scalar

def einstein_tensor(ricci_t, ricci_s, metric):
    pass

if __name__ == "__main__":
    from pprint import pprint
    
    t = sympy.Symbol('t')
    r = sympy.Symbol('r')
    theta = sympy.Symbol('theta')
    phi = sympy.Symbol('phi')
    k = sympy.Symbol('k')
    a = sympy.Function('a')(t)



    metric = np.diag([-1, a**2/(1-k*r**2), a**2*r**2,a**2*r**2*sympy.sin(theta)**2])
    metric_key = [t, r, theta, phi]
    c = christoffel(metric, metric_key)
    print("Christoffel symbols calculated")
    r = reimann_tensor(c, metric_key)
    print("Reimann tensor calculated")
    ri = ricci_tensor(r)
    s = ricci_scalar(ri, metric)
    print(s)
