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
import matplotlib.pyplot as pplot

def dVdt(V, t):
	return sqrt(3*I0*(omega['m']*V0*V+omega['r']*V0**sp.Rational(4, 3)*V**sp.Rational(2, 3)+omega['v']*V**2)+c)

def make_dHdt(V):
	def dHdt(H, t):
		V_t = get_value(V, t)
		return (I0/2*(omega['m']*V0+sp.Rational(2, 3)*omega['r']*V0**sp.Rational(4, 3)*V_t**sp.Rational(-1, 3)+2*omega['v']*V_t)-H*dVdt(V_t,t))/V_t
	return dHdt

def make_dSdt(H):
	def dSdt(S, t):
	    return get_value(H, t) * S
	return dSdt

def hubble_parameters():
	V = RK4(dVdt, V0, start, stop, step/4.0)
	dHdt = make_dHdt(V)
	return RK4(dHdt, A0, start, stop, step/2.0), RK4(dHdt, B0, start, stop, step/2.0), RK4(dHdt, C0, start, stop, step/2.0)

def scale_factors(Ha, Hb, Hc):
	dadt = make_dSdt(Ha)
	dbdt = make_dSdt(Hb)
	dcdt = make_dSdt(Hc)
	return RK4(dadt, a0, start, stop, step), RK4(dbdt, b0, start, stop, step), RK4(dcdt, c0, start, stop, step)

def plot_hubble_parameters():
	Ha, Hb, Hc = hubble_parameters()
	pplot.scatter(t, np.float64(values_at_times(Ha, t)), c = 'r')
	pplot.scatter(t, np.float64(values_at_times(Hb, t)), c = 'g')
	pplot.scatter(t, np.float64(values_at_times(Hc, t)), c = 'b')
	pplot.title('Hubble Parameters')
	pplot.show()
	return Ha, Hb, Hc

def plot_scale_factors(Ha, Hb, Hc):
	a, b, c = scale_factors(Ha, Hb, Hc)
	pplot.scatter(t, np.float64(values_at_times(a, t)), c = 'r')
	pplot.scatter(t, np.float64(values_at_times(b, t)), c = 'g')
	pplot.scatter(t, np.float64(values_at_times(c, t)), c = 'b')
	pplot.title('Scale Factors')
	pplot.show()

def RK4(dfdt, f0, start, stop, step):
	f = {start: f0}
	t = start
	val = f0
	cond = iter_cond(start, stop)
	while (cond(t, stop)):
		k1 = dfdt(val, t)
		k2 = dfdt(val + step/2.0*k1, t + step/2.0)
		k3 = dfdt(val + step/2.0*k2, t + step/2.0)
		k4 = dfdt(val + step*k3, t + step)
		t = t + step
		val = val + step/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
		set_value(f, t, val)
	return f

def iter_cond(start, stop):
	if (start < stop):
		def cond(a, b):
			return round(a, 8) < round(b, 8)
	else:
		def cond(a, b):
			return round(a, 8) > round(b, 8)
	return cond

def set_value(f, t, v):
	f[round(t, 8)] = v

def get_value(f, t):
	return f[round(t, 8)]

def values_at_times(v, t):
	values = []
	for time in t:
		values.append(get_value(v, time))
	return values

if __name__ == '__main__':
    #Initial directional Hubble constants
    A0 = 1.0
    B0 = 2.0
    C0 = 3.0

    #Initial directional scale factors
    a0 = 0
    b0 = 0
    c0 = 0

    #Farctional energy-densities of the universe
    omega = {'m': 0, 'r': 1, 'v': 0}

    #Times at which to calculate functions
    start = 0
    stop = 1
    step = 0.05

    I0 = A0*B0+A0*C0+B0*C0
    H0 = A0+B0+C0
    V0 = a0*b0*c0
    c = V0**2*(H0**2-3*I0)
    t = np.linspace(start, stop, abs((stop-start)/step)+1)

    Ha, Hb, Hc = plot_hubble_parameters()
    plot_scale_factors(Ha, Hb, Hc)
    