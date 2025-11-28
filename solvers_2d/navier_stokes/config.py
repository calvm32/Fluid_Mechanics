from firedrake import *

# constants
T = 2           # final time
dt = 0.1        # timestepping length
theta = 1/2     # theta constant
N = 10          # mesh size

Re = Constant(1.0)    # Reynold's num for viscosity

# functions
ufl_v = as_vector([1, 0])           # velocity ic
ufl_p = Constant(0.0)               # pressure ic
ufl_f = as_vector([0, 0])           # source term f
ufl_g = as_vector([0, 0])           # bdy condition g

# -------- 
# For MMS 
# --------

P = Constant(10.0)      # pressure constant
H = 5.0                 # height of rectangle, just take length = 3H

t = Constant(0.0) # symbolic constant for t
ufl_exp = ufl.exp # ufl e, so t gets calculated correctly

# exact functions for Poiseuille flow 
ufl_v_exact = as_vector(            # velocity ic
    [Re*( sin(pi*y/H)*
         ufl_exp(((pi**2)*t)/(H**2)) + 
         0.5*P*y**2 + 0.5*P*H*y ), 
     Constant(0.0)])
ufl_p_exact = P                     # pressure ic
ufl_f_exact = as_vector([1, 0])     # source term f
ufl_g_exact = as_vector([1, 0])     # bdy condition g