import scipy.integrate
import numpy
import matplotlib.pyplot

"""
Represents a universe described by three directional scale factors whose evolution is described by second order differential equations

This is an abstract class -- it needs a dydt implimentation to work

Subclasses must implement dydt, and should be sure to call the superclass constructor

Subclasses might also want to implement more variable calculations, which can be done by overriding calculate_variable
Just be sure to default to the superclass method if the subclass doesn't implement a variable
"""
class AnisotropicUniverse:

   def __init__(self):
      self.var_names = {'sf': 'Scale Factors', 'v': 'Volume', 'sfd': 'Scale Factor Derivatives', 'h': 'Hubble Parameters', 'hr': 'Hubble Parameter Ratios'}

   #Initial condiitons are taken to be specified at the time closest to 0 in times
   #dydt must return values in the form [a_dot, a_dot_dot, b_dot, b_dot_dot, c_dot, c_dot_dot]
   def evolve(self, initial_conditions, times):
      self.initial_conditions = initial_conditions
      self.a0, self.a_dot0, self.b0, self.b_dot0, self.c0, self.c_dot0 = initial_conditions
      self.V0 = self.a0*self.b0*self.c0
      self.times = times
      t0_index = min((abs(val), idx) for (idx, val) in enumerate(times))[1]
      if t0_index == 0:
         self.data = scipy.integrate.odeint(self.dydt, initial_conditions, times).transpose()
      else:
         first_half = scipy.integrate.odeint(self.dydt, initial_conditions, times[:t0_index+1][::-1])[::-1]
         second_half = scipy.integrate.odeint(self.dydt, initial_conditions, times[t0_index:])[1:]
         self.data = numpy.concatenate((first_half, second_half)).transpose()

   #var is an abbreviation for a varibale
   #Options are listed in the array in the constructor
   def calculate_variable(self, var):
      data = self.data
      if var == 'sf':
         return [data[0], data[2], data[4]]
      elif var == 'v':
         return [[data[0][i]*data[2][i]*data[4][i] for i in range(len(self.times))]]
      elif var == 'sfd':
         return [data[1], data[3], data[5]]
      elif var == 'h':
         return [[data[j][i]/data[j-1][i] for i in range(len(self.times))] for j in [1, 3, 5]]
      elif var == 'hr':
         Ha, Hb, Hc = self.calculate_variable('h')
         return [[H1[i]/H2[i] for i in range(len(self.times))] for [H1, H2] in [[Hb, Ha], [Ha, Hc], [Ha, Hb]]]

   #varaibles is a list of variable abbreviations
   #Options are listed in the array in the constructor
   def print_variables(self, variables):
      for var in variables:
         print(self.var_names[var])
         print(self.calculate_variable(var))
         print('\n')

   #varaibles is a list of variable abbreviations
   #Options are listed in the array in the constructor
   def plot_variables(self, variables):
      colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
      for var in variables:
         matplotlib.pyplot.title(self.var_names[var])
         variable = self.calculate_variable(var)
         c = 0
         for l in variable:
            matplotlib.pyplot.scatter(self.times, l, c = colors[c])
            c = c+1
         matplotlib.pyplot.show()





