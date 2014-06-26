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

#w = sp.Rational(1, 3) #Determines your fluid
w = sp.Rational(1, 3)

I = A*B+A*C+B*C
H = A+B+C

#c = .5#Integration constant, equal to sum of Ks
#V0 = sp.sqrt(c/(H**2*(1-(3*I/H**2)))) #Initial volume

V0 = 0.01 #Initial volume
c = V0**2*(H**2*(1-(3*I/H**2)))

def dVdt(v, t):
    return sqrt(3*I*V0**(1+w)*v**(1-w)+c)

times = np.linspace(0, 2, 100)
V = scipy.integrate.odeint(dVdt, V0, times)
Vdict = dict(zip(times, list(V.T[0])))

def dHdt(h, t):
    return (I/2*Vdict[t]**(-w)*V0**(1+w)*(1-w)-h*dVdt(Vdict[t],t))/Vdict[t]

def euler(dfdt, initial_condition, times):
    vals = [initial_condition]
    for n,t in enumerate(times[:-1]):
        vals.append(vals[-1]+dfdt(vals[-1], t)*(times[n+1]-t))
    return vals

Ha = np.float64(euler(dHdt, A, times))
print(Ha)

Hb = np.float64(euler(dHdt, B, times))
print(Hb)

Hc = np.float64(euler(dHdt, C, times))
print(Hc)

pplot.scatter(times, Hb, c = 'b')
pplot.scatter(times, Ha, c = 'r')
pplot.scatter(times, Hc, c = 'g')

#pplot.scatter(times, V)
#pplot.scatter(times, [((t+.5)/.5)**(3.0/2) for t in times])
#print(list(V.T[0])[-1])
pplot.show()