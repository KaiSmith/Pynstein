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
    symbols = [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]
    for alpha in range(4):
        for beta in range(4):
            for gamma in range(4):
                total = 0
                for delta in range(4):
                    total += 0.5*np.matrix(metric).I[alpha][delta] * (sympy.diff(metric[delta][beta], metric_key[gamma]) + 
                            sympy.diff(metric[delta][gamma], metric_key[beta]) + sympy.diff(metric[beta][gamma], metric_key[delta]))
                symbols[alpha][beta][gamma] = total
                    

def reimann_tensor(chris):
    pass

def ricci_tensor(reimann):
    pass

def ricci_scalar(ricci_t, metric):
    pass

def einstein_tensor(ricci_t, ricci_s, metric):
    pass

if __name__ == "__main__":
    t = sympy.Symbol('t')
    r = sympy.Symbol('r')
    theta = sympy.Symbol('theta')
    phi = sympy.Symbol('phi')
    k = sympy.Symbol('k')
    a = sympy.Function('a')(t)



    metric = np.array([[-1, 0, 0, 0],[0, a**2/(1-k*r**2), 0, 0],[0, 0, a**2*r**2, 0],[0, 0, 0, a**2*r**2*sympy.sin(theta)**2]])
    metric_key = [t, r, theta, phi]
    c = christoffel(metric, metric_key)
    from pprint import pprint
    pprint(c)
