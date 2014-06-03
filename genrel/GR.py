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
    for alpha in range(4):
        for beta in range(4):
            for gamma in range(4):
                total = 0
                for delta in range(4):                    
                    total += sympy.Matrix(metric).inv()[alpha, delta] * (sympy.diff(metric[delta][beta], metric_key[gamma]) + 
                            sympy.diff(metric[delta][gamma], metric_key[beta]) - sympy.diff(metric[beta][gamma], metric_key[delta]))
                symbols[alpha][beta][gamma] = sympy.cancel(1*total/2)
    return symbols
                    

def reimann_tensor(chris):
    pass

def ricci_tensor(reimann):
    pass

def ricci_scalar(ricci_t, metric):
    pass

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
    pprint(c)
