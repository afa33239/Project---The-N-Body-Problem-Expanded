## holds configurations of the simulation & simulation itself and runs it


import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from typing import List

from code.nbody.bodies import Body, SystemState

from code.nbody.integrators.euler import EulerIntegrator

from code.nbody.solvers.direct import DirectSolver

from code.nbody.physics import compute_kinetic_energy, compute_potential_energy, compute_angular_momentum


class SimulationConfig: 

    def __init__(self, dt, timesteps, softening=0.001):
        self.dt = dt
        self.timesteps = timesteps
        self.softening = softening


class Simulation:

    def __init__(self, bodies: List[Body], cfg: SimulationConfig, integrator=None, solver=None):
        self.state = SystemState(bodies)
        self.cfg = cfg
        self.integrator = integrator or EulerIntegrator()
        self.solver = solver or DirectSolver()

        self.energy_history = []  #keeps track of energies and angular momentum
        self.kinetic_history = []
        self.potential_history = []
        self.angular_momentum_history = []

        self.energy_drift = [] #this is used to track relative energy drift over time 
        self.angular_momentum_drift = []


    def run(self): #coordinates the simulation
        self._clear_histories()


        accel_fn = self._initialize_simulation()

        pss = []
        ps = [(b.x, b.y) for b in self.state.bodies]
        pss.append(ps)



        for _ in range(self.cfg.timesteps):
            self._step(accel_fn)
            ps = [(b.x, b.y) for b in self.state.bodies]
            pss.append(ps)

        return pss
    


    def _initialize_simulation(self):
        def accel_fn(bodies):
            return self.solver.accelerations(bodies, self.cfg)

        # convert (x0, v0) -> leapfrog-native (x0, v1/2)
        self.state = self.integrator.initialize(self.state, self.cfg, accel_fn)

        # diagnostic, time-aligned state (x0, v0) for measurements
        diag = self.integrator.synchronize(self.state, self.cfg, accel_fn)

        # store initial diagnostic snapshot
        self.state_history = [diag.copy()]

        # compute initial invariants from diagnostic state
        K0 = compute_kinetic_energy(diag.bodies)
        U0 = compute_potential_energy(diag.bodies, self.cfg)
        E0 = K0 + U0
        L0 = compute_angular_momentum(diag.bodies)

        self.E0 = E0
        self.L0 = L0

        self.kinetic_history.append(K0)
        self.potential_history.append(U0)
        self.energy_history.append(E0)
        self.energy_drift.append(0.0)

        self.angular_momentum_history.append(L0)
        self.angular_momentum_drift.append(0.0)

        return accel_fn
    


    def _step(self, accel_fn):
        # advance integrator internal state by one step
        self.state = self.integrator.step(self.state, self.cfg, accel_fn)

        # convert to diagnostic (time-aligned) state for history + invariants
        diag = self.integrator.synchronize(self.state, self.cfg, accel_fn)
        self.state_history.append(diag.copy())

        # compute invariants from diagnostic state
        K = compute_kinetic_energy(diag.bodies)
        U = compute_potential_energy(diag.bodies, self.cfg)
        E = K + U
        L = compute_angular_momentum(diag.bodies)

        self.kinetic_history.append(K)
        self.potential_history.append(U)
        self.energy_history.append(E)
        self.energy_drift.append((E - self.E0) / abs(self.E0))

        self.angular_momentum_history.append(L)
        self.angular_momentum_drift.append((L - self.L0) / abs(self.L0))


    def _clear_histories(self): #method to clean the histories before running the simulation
        self.state_history = []
        self.energy_history = []
        self.kinetic_history = []
        self.potential_history = []
        self.angular_momentum_history = []
        self.energy_drift = []
        self.angular_momentum_drift = [] 


    def closestDistance(self):
        import math

        ret = None
        pss = self.run()

        for step in pss:
            for a in step:
                for b in step:
                    if a is not b:
                        d = a.squareDist(b)
                        if ret is None or d < ret:
                            ret = d

        return None if ret is None else math.sqrt(ret)

    def show(self, x0, y0, x1, y1):
       
        pss = self.run()

        fig, ax = plt.subplots()
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)

        scatters = [
            ax.scatter([], [], label=f"Body {i}")
            for i in range(len(self.state.bodies))
        ]

        time_text = ax.text(
            0.02, 0.98, '',
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        )

        def update(frame):
            for i, sc in enumerate(scatters):
                sc.set_offsets([pss[frame][i].x, pss[frame][i].y])
            time_text.set_text(f"Timestep: {frame}")
            return scatters + [time_text]

        FuncAnimation(
            fig, update,
            frames=len(pss),
            interval=1,
            blit=True,
            repeat=False
        )

        ax.set_xlabel("X coordinate (AU)")
        ax.set_ylabel("Y coordinate (AU)")
        ax.set_title("Celestial Body Trajectories")
        ax.legend()
        plt.show()