from firedrake import *

def create_timestep_solver(get_data, dsN, theta, u_old, u_new, make_weak_form,
                           bcs, nullspace, solver_parameters, appctx, W):
    """
    Prepare timestep solver by theta-scheme for given
    function get_data(t) returning data (f(t), g(t)), given
    solution u_old at time t and unknown u_new at time t + dt.
    Return a solve function taking (t, dt).
    """

    # Default solver settings
    solver_kwargs = {}
    if bcs is not None:
        solver_kwargs["bcs"] = bcs
    if nullspace is not None:
        solver_kwargs["nullspace"] = nullspace
    if solver_parameters is not None:
        solver_kwargs["solver_parameters"] = solver_parameters
    if appctx is not None:
        solver_kwargs["appctx"] = appctx
    if W is not None:
        solver_kwargs["W"] = W

    # Initialize coefficients
    f_n, g_n = get_data(0)
    f_np1, g_np1 = get_data(0)
    idt = Constant(1.0)

    # Extract function space
    Z = u_new.function_space()

    # callable weak form
    weak_form = make_weak_form(theta, idt, f_n, f_np1, g_n, g_np1, dsN)

    def solve_(t, dt):
        """
        Update problem data to interval (t, t+dt) and run solver
        """

        # Update coefficients to current t, dt
        get_data(t, (f_n, g_n))
        get_data(t+dt, (f_np1, g_np1))
        idt.assign(1/dt)

        # split unknowns and previous solutions (symbolic split for form construction)
        u_sym, p_sym = split(u_new)    # symbolic unknowns (UFL)
        u_old_sym, p_old_sym = split(u_old)  # symbolic previous step (UFL)

        # build the residual using the weak form factory
        weak_form = make_weak_form(theta, idt, f_n, f_np1, g_n, g_np1, dsN)
        F = weak_form(u_sym, p_sym, u_old_sym, p_old_sym, *TestFunctions(u_new.function_space()))

        # Solve the nonlinear problem F==0 for u_new
        # Note: solve() will compute the Jacobian automatically (derivative w.r.t. u_new),
        # because the form only contains split(u_new) symbolic UFL objects, not subfunctions.
        solve(F == 0, u_new, **solver_kwargs)

    return solve_