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
from math import pi

a0 = 1.0
b0 = 1.0
c0 = 1.0

a_dot0 = 1.0
b_dot0 = 1.0
c_dot0 = 1.0

A0 = a_dot0/a0
B0 = b_dot0/b0
C0 = c_dot0/c0

omega0 = 1

#Open  -1
#Flat   0
#Closed 1
k = 1

t = np.linspace(0, 100, 1000)

I0 = A0*B0+B0*C0+A0*C0
H0 = A0+B0+C0
V0 = a0*b0*c0

chi0 = (omega0*I0*H0)/(3*(a_dot0+b_dot0+c_dot0))

#const = 8*pi*G*p0
const = 1

def dydt(y, t):
	a, a_dot, b, b_dot, c, c_dot = y
	"""
	a_dot_dot = (a/2.0)*(chi0*(a0/a - b0/b - c0/c)*(V0/(a*b*c) + k) - (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) + (b_dot*c_dot)/(b*c))
	b_dot_dot = (b/2.0)*(chi0*(-a0/a + b0/b - c0/c)*(V0/(a*b*c) + k) - (a_dot*b_dot)/(a*b) + (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
	c_dot_dot = (c/2.0)*(chi0*(-a0/a - b0/b + c0/c)*(V0/(a*b*c) + k) + (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
	"""
	a_dot_dot = (a/2.0)*(-const*(V0/(a*b*c))*(-a0/a + b0/b + c0/c) - k*(-1/a**2 + 1/b**2 + 1/c**2) - (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) + (b_dot*c_dot)/(b*c))
	b_dot_dot = (b/2.0)*(-const*(V0/(a*b*c))*(a0/a - b0/b + c0/c) -k*(1/a**2 - 1/b**2 + 1/c**2) - (a_dot*b_dot)/(a*b) + (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
	c_dot_dot = (c/2.0)*(-const*(V0/(a*b*c))*(a0/a + b0/b - c0/c) -k*(1/a**2 + 1/b**2 - 1/c**2) + (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
	return [a_dot, a_dot_dot, b_dot, b_dot_dot, c_dot, c_dot_dot]

def plot_evolution():
	t = np.linspace(1, 5, 100)
	y0 = [a0, a_dot0, b0, b_dot0, c0, c_dot0]
	y = scipy.integrate.odeint(dydt, y0, t)

	a = [value[0] for value in y]
	a_dot = [value[1] for value in y]
	b = [value[2] for value in y]
	b_dot = [value[3] for value in y]
	c = [value[4] for value in y]
	c_dot = [value[5] for value in y]

	stop = len(t) - 1
	for values in [a, a_dot, b, b_dot, c, c_dot]:
		for i in range(1, len(t)):
			if abs(values[i]/values[i-1]) > 1000 and i < stop:
				stop = i
				break
	a, a_dot, b, b_dot, c, c_dot, t = a[:stop], a_dot[:stop], b[:stop], b_dot[:stop], c[:stop], c_dot[:stop], t[:stop]

	A = [a_dot[i]/a[i] for i in range(len(t))]
	B = [b_dot[i]/b[i] for i in range(len(t))]
	C = [c_dot[i]/c[i] for i in range(len(t))]

	V = [a[i]*b[i]*c[i] for i in range(len(t))]

	"""
	pplot.scatter(t, a_dot, c = 'r')
	pplot.scatter(t, b_dot, c = 'g')
	pplot.scatter(t, c_dot, c = 'b')
	pplot.title('First Derivatives')
	pplot.show()
	"""

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

	pplot.scatter(t, V, c = 'r')
	pplot.title('Volume')
	pplot.show()


def print_long_term_ratios():
	t = np.linspace(0, 1000000, 100000)
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

plot_evolution()