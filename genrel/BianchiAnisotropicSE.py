"""
Numerically evolves a Bianchi Class I universe with anisotripic stress energy given initial conditions
David Clark, Kai Smith
Case Western Reserve University
2014
"""
from math import *
import numpy as np
import sympy as sp
import scipy.integrate
import matplotlib.pyplot as pplot

a0 = 0.5
b0 = 1.0
c0 = 1.5

a_dot0 = 1.0
b_dot0 = 1.5
c_dot0 = 2.0

p0 = 1

p1 = a0**2*b0*c0*p0
p2 = a0*b0**2*c0*p0
p3 = a0*b0*c0**2*p0

t = np.linspace(0, 10, 50)

def dydt(y, t):
	a, a_dot, b, b_dot, c, c_dot = y
	a_dot_dot = ((a*b_dot*c_dot)/(b*c) - (a_dot*c_dot)/c - (a_dot*b_dot)/b + p3/(b*c**2)+ p2/(b**2*c) - p1/(a*b*c))/2.0
	b_dot_dot = ((b*a_dot*c_dot)/(a*c) - (b_dot*a_dot)/a - (b_dot*c_dot)/c + p3/(a*c**2) + p1/(a**2*c) - p2/(a*b*c))/2.0
	c_dot_dot = ((c*a_dot*b_dot)/(a*b) - (c_dot*a_dot)/a - (c_dot*b_dot)/b + p1/(a**2*b) + p2/(a*b**2) - p3/(a*b*c))/2.0
	return [a_dot, a_dot_dot, b_dot, b_dot_dot, c_dot, c_dot_dot]


def plot_stuff():
	y0 = [a0, a_dot0, b0, b_dot0, c0, c_dot0]
	y = scipy.integrate.odeint(dydt, y0, t)

	a = [value[0] for value in y]
	b = [value[2] for value in y]
	c = [value[4] for value in y]

	A = [value[1]/value[0] for value in y]
	B = [value[3]/value[2] for value in y]
	C = [value[5]/value[4] for value in y]

	V = [value[0]*value[2]*value[4] for value in y]
	#V_dot = [value[1]*value[2]*value[4] + value[0]*value[3]*value[4] + value[0]*value[2]*value[5] for value in y]

	pplot.scatter(t, np.float64(A), c = 'r')
	pplot.scatter(t, np.float64(B), c = 'g')
	pplot.scatter(t, np.float64(C), c = 'b')
	pplot.show()

	pplot.scatter(t, np.float64(a), c = 'r')
	pplot.scatter(t, np.float64(b), c = 'g')
	pplot.scatter(t, np.float64(c), c = 'b')
	pplot.show()

plot_stuff()