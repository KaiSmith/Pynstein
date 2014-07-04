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
	x1, x2, x3, x4, x5, x6 = y
	x1_dot = x2
	x2_dot = ((x1*x4*x6)/(x3*x5) - (x2*x6)/x5 - (x2*x4)/x3 + p3/(x3*x5**2)+ p2/(x3**2*x5) - p1/(x1*x3*x5))/2.0
	x3_dot = x4
	x4_dot = ((x3*x2*x6)/(x1*x5) - (x4*x2)/x1 - (x4*x6)/x5 + p3/(x1*x5**2) + p1/(x1**2*x5) - p2/(x1*x3*x5))/2.0
	x5_dot = x6
	x6_dot = ((x5*x2*x4)/(x1*x3) - (x6*x2)/x1 - (x6*x4)/x3 + p1/(x1**2*x3) + p2/(x1*x3**2) - p3/(x1*x3*x5))/2.0
	return [x1_dot, x2_dot, x3_dot, x4_dot, x5_dot, x6_dot]


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
	V_dot = [value[1]*value[2]*value[4] + value[0]*value[3]*value[4] + value[0]*value[2]*value[5] for value in y]

	pplot.scatter(t, np.float64(A), c = 'r')
	pplot.scatter(t, np.float64(B), c = 'g')
	pplot.scatter(t, np.float64(C), c = 'b')
	pplot.show()

	pplot.scatter(t, np.float64(a), c = 'r')
	pplot.scatter(t, np.float64(b), c = 'g')
	pplot.scatter(t, np.float64(c), c = 'b')
	pplot.show()

plot_stuff()