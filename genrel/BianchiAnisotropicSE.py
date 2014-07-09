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
b_dot0 = 1.0
c_dot0 = 1.0

A0 = a_dot0/a0
B0 = b_dot0/b0
C0 = c_dot0/c0

omega0 = 1

t = np.linspace(0, 1000000000, 1000)

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

	"""
	g_prime = [p0*(a0/value[0] - b0/value[2]) for value in y]
	f_prime = [rho0*((a0*b0*c0)/(value[0]*value[2]*value[4]))**(1.0/3.0) - p0*(c0/value[4]) for value in y]
	ratio = [g_prime[i]/f_prime[i] for i in range(0, len(g_prime))]
	limit = [(1 - ratio[i])/(1 + ratio[i]) for i in range(len(ratio))]
	print(limit[-1])
	"""

	"""
	A_over_B = [(value[1]/value[0])/(value[3]/value[2]) for value in y]
	print(A_over_B)
	"""

	"""
	step_size = t[1] - t[0]
	f = [0]
	for i in range(1, len(f_prime)):
		f.append(f[i-1] + f_prime[i]*step_size)
	g = [0]
	for i in range(1, len(g_prime)):
		g.append(g[i-1] + g_prime[i]*step_size)
	num = [f[i] + g[i] for i in range(len(f))]
	denom = [f[i] - g[i] for i in range(len(f))]
	ratio = [num[i]/denom[i] for i in range(1, len(f))]

	pplot.scatter(t[1:], np.float64(ratio), c = 'r')
	pplot.show()
	"""

	"""
	g_over_f = [g[i]/f[i] for i in range(1, len(f))]

	pplot.scatter(t[1:], np.float64(g_over_f), c = 'r')
	pplot.show()
	"""

	"""
	B_over_C = [A[i]/B[i] for i in range(len(t))]
	C_over_A = [C[i]/A[i] for i in range(len(t))]
	B_over_A = [B[i]/A[i] for i in range(len(t))]
	
	pplot.scatter(t, A, c = 'r')
	pplot.scatter(t, B, c = 'g')
	pplot.scatter(t, C, c = 'b')
	pplot.show()

	pplot.scatter(t, a, c = 'r')
	pplot.scatter(t, b, c = 'g')
	pplot.scatter(t, c, c = 'b')
	pplot.show()
	"""

	gamma0 = H0/(3.0*(a_dot0+b_dot0+c_dot0))
	f_prime = [(V0/V[i])**(1.0/3.0) - gamma0*(c0/c[i]) for i in range(len(t))]
	g_prime = [gamma0*(a0/a[i] - b0/b[i]) for i in range(len(t))]
	ratio = [g_prime[i] - f_prime[i] for i in range(len(t))]
	print(ratio[-1])

	"""
	print(frac)
	pplot.scatter(t, np.float64(frac), c = 'r')
	pplot.show()
	"""

plot_stuff()
