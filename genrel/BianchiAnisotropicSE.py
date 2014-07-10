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

a0 = 1.0
b0 = 1.0
c0 = 1.0

a_dot0 = 1.0
b_dot0 = 5.0
c_dot0 = 100000.0

A0 = a_dot0/a0
B0 = b_dot0/b0
C0 = c_dot0/c0

omega0 = 1

t = np.linspace(0, 2, 50)

I0 = A0*B0+B0*C0+A0*C0
H0 = A0+B0+C0
V0 = a0*b0*c0

chi0 = (omega0*I0*V0*H0)/(3*(a_dot0+b_dot0+c_dot0))

def dydt(y, t):
	a, a_dot, b, b_dot, c, c_dot = y
	a_dot_dot = (a/2.0)*((chi0/(a*b*c))*(a0/a - b0/b - c0/c) - (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) + (b_dot*c_dot)/(b*c))
	b_dot_dot = (b/2.0)*((chi0/(a*b*c))*(-a0/a + b0/b - c0/c) - (a_dot*b_dot)/(a*b) + (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
	c_dot_dot = (c/2.0)*((chi0/(a*b*c))*(-a0/a - b0/b + c0/c) + (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
	return [a_dot, a_dot_dot, b_dot, b_dot_dot, c_dot, c_dot_dot]

def plot_evolution():
	y0 = [a0, a_dot0, b0, b_dot0, c0, c_dot0]
	y = scipy.integrate.odeint(dydt, y0, t)

	a = [value[0] for value in y]
	b = [value[2] for value in y]
	c = [value[4] for value in y]

	A = [value[1]/value[0] for value in y]
	B = [value[3]/value[2] for value in y]
	C = [value[5]/value[4] for value in y]

	B_over_C = [A[i]/B[i] for i in range(len(t))]
	C_over_A = [C[i]/A[i] for i in range(len(t))]
	B_over_A = [B[i]/A[i] for i in range(len(t))]
	
	pplot.scatter(t, a, c = 'r')
	pplot.scatter(t, b, c = 'g')
	pplot.scatter(t, c, c = 'b')
	pplot.title('Scale Factors')
	pplot.show()

	pplot.scatter(t, A, c = 'r')
	pplot.scatter(t, B, c = 'g')
	pplot.scatter(t, C, c = 'b')
	pplot.title('Hubble Parameters')
	pplot.show()

	pplot.scatter(t, B_over_C, c = 'r')
	pplot.scatter(t, C_over_A, c = 'g')
	pplot.scatter(t, B_over_A, c = 'b')
	pplot.title('Hubble Parameter Ratios')
	pplot.show()

def print_long_term_ratios():
	t = np.linspace(0, 1000000000, 10000000)
	y0 = [a0, a_dot0, b0, b_dot0, c0, c_dot0]
	y = scipy.integrate.odeint(dydt, y0, t)

	A = [value[1]/value[0] for value in y]
	B = [value[3]/value[2] for value in y]
	C = [value[5]/value[4] for value in y]

	B_over_C = [A[i]/B[i] for i in range(len(t))]
	C_over_A = [C[i]/A[i] for i in range(len(t))]
	B_over_A = [B[i]/A[i] for i in range(len(t))]

	print('B/C: ' + str(B_over_C[-1]))
	print('C/A: ' + str(C_over_A[-1]))
	print('B/A: ' + str(B_over_A[-1]))

#plot_evolution()
print_long_term_ratios()