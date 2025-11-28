from firedrake import *

# constants
T = 2           # final time
dt = 0.1        # timestepping length
theta = 1/2     # theta constant
N = 10          # mesh size

# functions
ufl_f = cos(x*pi)*cos(y*pi)     # source term f
ufl_g = 0                       # bdy condition g
ufl_u0 = 0                      # initial condition u0