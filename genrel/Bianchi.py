from GR import *
import numpy as np
import sympy as sp

import scipy.integrate
import matplotlib.pyplot as pplot
import pylab
from math import *

#Initial directional Hubble constants
A = 1.0 #0.98
B = 1.1 #1.00
C = -0.01 #1.02

w = {'m': 0, 'r': sp.Rational(1, 3), 'v': -1}
omega = {'m': 0.5, 'r': 0.5, 'v': 0}

I = A*B+A*C+B*C
H = A+B+C

V0 = 0.01 #Initial volume
c = V0**2*(H**2*(1-(3*I/H**2)))

def dVdt(v, t):
    return sqrt(3*I*(omega['m']*V0*v+omega['r']*V0**sp.Rational(4, 3)*v**sp.Rational(2, 3)+omega['v']*v**2)+c)

times = np.linspace(0, 2, 100)
V = scipy.integrate.odeint(dVdt, V0, times)
Vdict = dict(zip(times, list(V.T[0])))

def dHdt(h, t):
    return (I/2*(omega['m']*V0+sp.Rational(2, 3)*omega['r']*V0**sp.Rational(4, 3)*Vdict[t]**sp.Rational(-1, 3)+2*omega['v']*Vdict[t])
            -h*dVdt(Vdict[t],t))/Vdict[t]

def euler(dfdt, initial_condition, times):
    vals = [initial_condition]
    for n,t in enumerate(times[:-1]):
        vals.append(vals[-1]+dfdt(vals[-1], t)*(times[n+1]-t))
    return vals

Ha = np.float64(euler(dHdt, A, times))
#print(Ha)

Hb = np.float64(euler(dHdt, B, times))
#print(Hb)

Hc = np.float64(euler(dHdt, C, times))
#print(Hc)

pplot.scatter(times, Hb, c = 'b')
pplot.scatter(times, Ha, c = 'r')
pplot.scatter(times, Hc, c = 'g')

#pplot.scatter(times, V)
#pplot.scatter(times, [((t+.5)/.5)**(3.0/2) for t in times])
#print(list(V.T[0])[-1])
pplot.show()
