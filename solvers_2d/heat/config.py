from firedrake import *

# ---------
# Constants
# ---------

T = 2           # final time
dt = 0.1        # timestepping length
theta = 1/2     # theta constant
N = 10          # mesh size

# ----------------
# For single solve 
# ----------------

N = 10          # mesh size

# -------------
# For MMS solve
# -------------

# MMS loops over mesh sizes in this list
N_list = []
for exp in range(1, 10):
    N = 2**exp
    N_list.append(N)