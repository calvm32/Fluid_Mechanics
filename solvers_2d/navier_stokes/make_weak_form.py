from firedrake import *

from .config_constants import Re

def make_weak_form(theta, idt, f, f_old, g, g_old, dx, dsN):
    def residual(u_new, u_old, v):
        u, p = split(u_new)
        u_old_, p_old = split(u_old)
        v_u, v_p = split(v)

        # theta-scheme time derivative
        F = idt * inner(u - u_old_, v_u) * dx

        # convection-diffusion term
        F += theta * (1/Re * inner(grad(u), grad(v_u)) + inner(dot(grad(u), u), v_u) - p * div(v_u)) * dx
        F += (1 - theta) * (1/Re * inner(grad(u_old_), grad(v_u)) + inner(dot(grad(u_old_), u_old_), v_u) - p_old * div(v_u)) * dx

        # incompressibility
        F += div(u) * v_p * dx

        return F
    return residual