import numpy
from AnisotropicUniverse import AnisotropicUniverse
from CurvedUniverse import CurvedUniverse

class S2E1Universe(AnisotropicUniverse, CurvedUniverse):
	def __init__(self, shape_or_k = 0):
		AnisotropicUniverse.__init__(self)
		CurvedUniverse.__init__(self, shape_or_k)

	def evolve(self, initial_conditions, times):
		a0, a_dot0, b0, b_dot0 = initial_conditions[:-2]
		Ha0, Hb0 = a_dot0/a0, b_dot0/b0
		self.const = -(self.k/a0**2.0 + 2.0*Ha0*Hb0 + Ha0**2.0)
		AnisotropicUniverse.evolve(self, initial_conditions, times)

	def dydt(self, y, t):
		a, a_dot, b, b_dot = y[:-2]
		k = self.k
		a0, a_dot0, b0, b_dot0 = self.initial_conditions[:-2]
		V = a**2.0*b
		V0 = a0**2.0*b0
		const = self.const

		a_dot_dot = -(a/2.0)*((a_dot**2.0)/(a**2.0) + k/a**2.0 + const*(V0/V)*(b0/b))
		b_dot_dot = -(b/2.0)*(-(a_dot**2.0)/(a**2.0) - k/a**2.0 + 2.0*(a_dot*b_dot)/(a*b) + const*(V0/V)*(2*(a0/a) - (b0/b)))
		c_dot_dot = 0

		return [a_dot, a_dot_dot, b_dot, b_dot_dot, 0, 0]

if __name__ == "__main__":
   universe = S2E1Universe('f')
   universe.evolve([1.0, 1.0, 1.0, 1.0, 0, 0], numpy.linspace(0, 10, 100))
   universe.plot_variables(['sf'])





