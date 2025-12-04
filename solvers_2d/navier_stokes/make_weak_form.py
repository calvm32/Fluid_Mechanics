from firedrake import *

from .config_constants import Re

def make_weak_form(theta, idt, f, f_old, g, g_old, dx , dsN):
    """
    Returns func F(u, u_old, p, q, v), 
    which builds weak form
    using external coefficients
    """

    f_mid = theta * f + (1-theta) * f_old
    g_mid = theta * g + (1-theta) * g_old

    def F(u_new, u_old, v):
        u_new, p_new = split(u_new)
        u_old_, p_old = split(u_old)
        v_u, v_p = split(v)

        u_mid = theta * u_new + (1 - theta) * u_old_
    
        return (
            # Time derivative
            idt*inner(u_new - u_old_, v_u)*dx

            # Diffusion
            + (1/Re)*inner(grad(u_mid), grad(v_u))*dx

            # Convection — Crank–Nicolson (implicit midpoint)
            + inner(dot(u_mid, grad(u_mid)), v_u)*dx

            # Pressure
            - inner(p_new, div(v_u))*dx
            + inner(div(u_mid), v_p)*dx

            # Source, boundary
            - inner(f_mid, v_u)*dx
            - inner(g_mid, v_u)*dsN
        )
    
    return F