"""
Numerically evolves a Bianchi Class I universe given initial conditions.
David Clark, Kai Smith
Case Western Reserve University
2014
"""
from GR import *
from math import *
import numpy as np
import sympy as sp
import scipy.integrate
import matplotlib.pyplot as pplot

#Initial directional Hubble constants
A0 = 1.2
B0 = 1.0
C0 = 0.8

#Initial directional scale factors
a0 = 1.0
b0 = 1.0
c0 = 1.0

#Farctional energy-densities of the universe
omega = {'m': 1.0/3.0, 'r': 1.0/3.0, 'v': 1.0/3.0}

#Times at which to calculate functions
t = np.linspace(0, 3, 100)

I0 = A0*B0+A0*C0+B0*C0
H0 = A0+B0+C0
V0 = a0*b0*c0
c = V0**2*(H0**2-3*I0)

def dVdt(V, t0):
    return sqrt(3*I0*(omega['m']*V0*V+omega['r']*V0**sp.Rational(4, 3)*V**sp.Rational(2, 3)+omega['v']*V**2)+c)

def make_dHdt(V):
	def dHdt(H, t0):
	    return (I0/2*(omega['m']*V0+sp.Rational(2, 3)*omega['r']*V0**sp.Rational(4, 3)*V[t0]**sp.Rational(-1, 3)+2*omega['v']*V[t0])
	    	-H*dVdt(V[t0],t0))/V[t0]
	return dHdt

def make_dSdt(H):
	def dSdt(S, t0):
	    return H[t0] * S
	return dSdt

def hubble_parameters():
	V = dict(zip(t, scipy.integrate.odeint(dVdt, V0, t)))
	dHdt = make_dHdt(V)
	return euler(dHdt, A0, t), euler(dHdt, B0, t), euler(dHdt, C0, t)

def scale_factors():
	Ha, Hb, Hc = hubble_parameters()
	Ha = dict(zip(t, Ha))
	Hb = dict(zip(t, Hb))
	Hc = dict(zip(t, Hc))
	dadt = make_dSdt(Ha)
	dbdt = make_dSdt(Hb)
	dcdt = make_dSdt(Hc)
	return euler(dadt, a0, t), euler(dbdt, b0, t), euler(dcdt, c0, t)

def plot_hubble_parameters():
	Ha, Hb, Hc = hubble_parameters()
	pplot.scatter(t, np.float64(Ha), c = 'r')
	pplot.scatter(t, np.float64(Hb), c = 'g')
	pplot.scatter(t, np.float64(Hc), c = 'b')
	pplot.title('Hubble Parameters')
	pplot.show()

def plot_scale_factors():
	a, b, c = scale_factors()
	pplot.scatter(t, np.float64(a), c = 'r')
	pplot.scatter(t, np.float64(b), c = 'g')
	pplot.scatter(t, np.float64(c), c = 'b')
	pplot.title('Scale Factors')
	pplot.show()

def euler(dfdt, f0, t):
    vals = [f0]
    for n,t0 in enumerate(t[:-1]):
        vals.append(vals[-1]+dfdt(vals[-1], t0)*(t[n+1]-t0))
    return vals

plot_hubble_parameters()
plot_scale_factors()
