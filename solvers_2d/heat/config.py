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

# -------- 
# For MMS 
# --------

t = Constant(0.0) # symbolic constant for t
ufl_exp = ufl.exp # ufl e, so t gets calculated correctly

# exact functions for u=e^t*sin(pix)*cos(piy)
ufl_u_exact = ufl_exp(t)*cos(pi*x)*cos(pi*y)                # source term f
ufl_f_exact = (1+2*pi**2)*ufl_exp(t)*cos(pi*x)*cos(pi*y)    # bdy condition g
ufl_g_exact = 0                                             # initial condition u0