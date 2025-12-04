from firedrake import * 

from solvers_2d.timestepper import timestepper
from .make_weak_form import make_weak_form
from solvers_2d.printoff import blue

from .config_constants import t0, T, dt, theta, N, solver_parameters, appctx, vtkfile_name

mesh = UnitSquareMesh(N, N)
x, y = SpatialCoordinate(mesh)

dx = Measure("dx", domain=mesh)
ds = Measure("ds", domain=mesh)

V = VectorFunctionSpace(mesh, "CG", 2)
W = FunctionSpace(mesh, "CG", 1)
Z = V * W

bcs = [DirichletBC(Z.sub(0), Constant((1, 0)), (4,)),
       DirichletBC(Z.sub(0), Constant((0, 0)), (1, 2, 3))]

nullspace = MixedVectorSpaceBasis(
    Z, [Z.sub(0), VectorSpaceBasis(constant=True)])

def get_data(t):
    
    # functions
    ufl_v0 = as_vector([sin(pi*x)*t, cos(pi*y)])    # velocity ic
    ufl_p0 = Constant(5.0)                          # pressure ic
    ufl_f = as_vector([0, 0])                       # source term f
    ufl_g = as_vector([0, 0])                       # bdy condition g

    ufl_u0 = [ufl_v0, ufl_p0]

    # returns
    return {"ufl_u0": ufl_u0,
            "ufl_f": ufl_f,
            "ufl_g": ufl_g,}

timestepper(get_data, theta, 
            Z, dx, ds(1), 
            t0, T, dt,
            make_weak_form=make_weak_form,
            bcs=bcs, nullspace=nullspace,
            olver_parameters=solver_parameters,
            appctx=appctx, vtkfile_name=vtkfile_name)
