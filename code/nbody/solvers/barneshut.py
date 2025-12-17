from code.nbody.solvers import Solver
from code.nbody.physics import compute_accelerations

class BarnesHutSolver(Solver):
    def accelerations(self, bodies, cfg):
        return compute_accelerations(bodies, cfg) #uses physics module function to computer accelerations, then returns them
    