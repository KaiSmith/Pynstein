"""
genrel package for GR calculations
David Clark, Kai Smith, David Cyncynates
Case Western Reserve university
2014
"""

import numpy as np
import sympy

def christoffel(metric):
    pass

def reimann_tensor(chris):
    pass

def ricci_tensor(reimann):
    pass

def ricci_scalar(ricci_t, metric):
    pass

def einstein_tensor(ricci_t, ricci_s, metric):
    pass

if __name__ == "__main__":
    metric = np([[-1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]])
    c = christoffel(metric)
