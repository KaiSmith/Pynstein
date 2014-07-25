class CurvedUniverse():
	def __init__(self, shape_or_k = 0):
		self.shapes = {'o': -1, 'f': 0, 'c': 1}
		self.set_shape_or_k(shape_or_k)

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