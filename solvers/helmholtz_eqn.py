from firedrake import *

def helmholtz_eqn_solver(f, g, nu, V):

    """
    Solver for the indefinite Helmholtz equation used in wave problems:
        - nu*(div)^2 u + u = f      on Omega
        - (div u) n = g             on bdy(Omega)
    """

    # trial, test funcs
    u = TrialFunction(V)
    v = TestFunction(V)

    # RHS weak form, LHS weak form
    a = ( (-nu)*inner(grad(u), grad(v)) + inner(u, v)) * dx
    L = inner(f, v) * dx

    # Dirichlet BC
    bc = DirichletBC(V, g, "on_boundary")

    # solve
    u = Function(V)
    solve(a == L, u, bcs = bc, solver_parameters={'ksp_type': 'cg', 'pc_type': 'none'})

    VTKFile("helmholtz.pvd").write(u)