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
    return sympy.cancel(scalar)

def einstein_tensor(ricci_t, ricci_s, metric):
    einstein = np.empty((4, 4), dtype = type(sympy.Symbol('Test')*1))
    for alpha in range(4):
        for beta in range(4):
            einstein[alpha][beta] = sympy.cancel(ricci_t[alpha][beta] - 0.5*metric[alpha][beta]*ricci_s)
    return einstein

def readable_print(tensor, index = []):
    for n, entry in enumerate(tensor):
        if type(entry) != type(np.array([])):
            if entry != 0:
                print(str(index + [n])+" : ")
                sympy.pprint(entry)
        else:
            readable_print(entry, index + [n])

def raise_index(tensor, metric):
    inverse = np.array(sympy.Matrix(metric).inv())
    raised_form = np.dot(inverse, tensor)
    return raised_form

def simplify_tensor(tensor):
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
    chris = christoffel(metric, metric_key)
    reimann = reimann_tensor(chris, metric_key)
    ricci_t = ricci_tensor(reimann)
    ricci_s = ricci_scalar(ricci_t, metric)
    einstein = einstein_tensor(ricci_t, ricci_s, metric)
    readable_print(raise_index(einstein, metric))
