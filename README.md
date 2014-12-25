![alt tag](https://raw.githubusercontent.com/KaiSmith/Pynstein/master/pynstein_logo.png)

##What is Pynstein?
Pynstein is a python library that allows the user to easily do general relativity calculations.

##Features

* Given a metric, Pynstein can calculate:
  * Inerse metrics
  * Christoffel Symbols
  * Reimann Curvature Tensor
  * Ricci Curvature Tensor
  * Conservation equation
  * Einstein Tensor
  * and given the stress energy tensor, the Einstein Equations
* Can numerically evolve a universe given inital conditions

##Dependencies
* numpy
* scipy
* sympy
* matplotlib

##Example Usage
```python
#Define variables
t = sp.Symbol('t')
r = sp.Symbol('r')
theta = sp.Symbol('theta')
phi = sp.Symbol('phi')

#Define rho and p functions
w = sp.Symbol('w')
rho = sp.Function('rho')(t)
p = w*rho

#Create FRW metric
frw_metric = np.diag([-1, a**2/(1-k*r**2), a**2*r**2,a**2*r**2*sp.sin(theta)**2])
 frw_metric_key = [t, r, theta, phi]

#Generate Einstein tensor
einstein = einstein_tensor_from_scratch(frw_metric, frw_metric_key)
einstein = raise_one_index(einstein, frw_metric)

#Generate the 2 distinct Einstein field equations
ein_eq = einstein_equations(einstein, np.diag([-rho, p, p, p]))
```

##Authors
Kai Smith (kai.smith@case.edu)

David Clark (davidclark@berkeley.edu)
