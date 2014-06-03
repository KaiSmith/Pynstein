"""
genrel package for GR calculations
David Clark, Kai Smith, David Cyncynates
Case Western Reserve university
2014
"""

import numpy, sympy

class EinsteinTensor:

    def __init__(self, obj):
        self.obj = obj

    def __str__(self):
        return "Hello!"