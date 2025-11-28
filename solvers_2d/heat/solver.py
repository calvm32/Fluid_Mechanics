from firedrake import *
from solvers_2d.timestepper import timestepper
from .make_weak_form import make_weak_form
from .config import T, dt, theta, N, ufl_f, ufl_g, ufl_u0

# mesh
mesh = UnitSquareMesh(N, N)
x, y = SpatialCoordinate(mesh)

# declare function space and interpolate functions
V = FunctionSpace(mesh, "CG", 1)

f = Function(V)
g = Function(V)
u0 = Function(V)

u0.interpolate(ufl_u0)

# make data for iterative time stepping
def get_data(t, result=None):
    """Create or update data"""
    if result is None: # only allocate memory if hasn't been yet
        f = Function(V)
        g = Function(V)
    else:
        f, g = result

    f.interpolate(ufl_f)
    g.interpolate(ufl_g)
    return f, g

# run
timestepper(V, ds(1), theta, T, dt, u0, get_data, make_weak_form)