

class Integrator:
    def initialize(self, state, cfg, accel_fn):
        # default: do nothing
        return state

    def step(self, state, cfg, accel_fn):
        raise NotImplementedError()

    def synchronize(self, state, cfg, accel_fn):
        """
        Return a state suitable for diagnostics where positions and velocities
        correspond to the same time level (x^n, v^n).
        Default: state already synchronized.
        """
        return state