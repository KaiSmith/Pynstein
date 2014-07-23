import numpy
from AnisotropicUniverse import AnisotropicUniverse

"""
Represents a radiation-filled Bianchi Type-1 universe with anisotripic stress-energy and optional curvature
"""

class BianchiRadiationUniverse(AnisotropicUniverse):
	def __init__(self, shape_or_k = 0, const = 1.0):
		AnisotropicUniverse.__init__(self)
		self.shapes = {'o': 1, 'f': 0, 'c': -1}
		self.set_shape_or_k(shape_or_k)
		self.const = const

	def set_shape_or_k(self, shape_or_k):
		if type(shape_or_k) == str:
			self.__dict__['shape'] = shape_or_k
			self.__dict__['k'] = self.shapes[shape_or_k]
		else:
			self.__dict__['k'] = shape_or_k
			self.__dict__['shape'] = [k for shape, k in self.shapes.items() if k == shape_or_k][0]

	def __setattr__(self, name, value):
		if name == 'shape' or name == 'k':
			self.set_shape_or_k(value)
		else:
			self.__dict__[name] = value


	def dydt(self, y, t):
		a, a_dot, b, b_dot, c, c_dot = y
		const = self.const
		k = self.k
		a0, a_dot0, b0, b_dot0, c0, c_dot0 = self.initial_conditions
		V0 = a0*b0*c0
		a_dot_dot = (a/2.0)*(-const*(V0/(a*b*c))*(-a0/a + b0/b + c0/c) - k*(-1/a**2 + 1/b**2 + 1/c**2) - (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) + (b_dot*c_dot)/(b*c))
		b_dot_dot = (b/2.0)*(-const*(V0/(a*b*c))*(a0/a - b0/b + c0/c) -k*(1/a**2 - 1/b**2 + 1/c**2) - (a_dot*b_dot)/(a*b) + (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
		c_dot_dot = (c/2.0)*(-const*(V0/(a*b*c))*(a0/a + b0/b - c0/c) -k*(1/a**2 + 1/b**2 - 1/c**2) + (a_dot*b_dot)/(a*b) - (a_dot*c_dot)/(a*c) - (b_dot*c_dot)/(b*c))
		return [a_dot, a_dot_dot, b_dot, b_dot_dot, c_dot, c_dot_dot]

if __name__ == "__main__":
   universe = BianchiRadiationUniverse()
   universe.evolve([1.0, 1.0, 1.0, 2.0, 1.0, 3.0], numpy.linspace(-0.1, 5, 100))
   universe.plot_variables(['sf', 'h'])





