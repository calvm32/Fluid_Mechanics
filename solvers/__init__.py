"""
solvers package: provides time-stepping for PDEs

Modules:
    - timestepper.py: fixed-step theta-scheme time integrator
    - timestepper_adaptive.py: adaptive time-stepping version
"""

from .timestepper import timestepper
from .timestepper_adaptive import timestepper_adaptive

__all__ = [
    "timestepper",
    "timestepper_adaptive",
]