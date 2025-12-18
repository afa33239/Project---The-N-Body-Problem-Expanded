import time
import random
import matplotlib.pyplot as plt

from code.nbody.bodies import Body
from code.nbody.engine import Simulation, SimulationConfig
from code.nbody.integrators.leapfrog import LeapfrogIntegrator
from code.nbody.solvers.direct import DirectSolver
from code.nbody.solvers.barneshut import BarnesHutSolver


# Utility functions

def make_random_bodies(
    N,
    seed=42,
    pos_scale=1.0,
    vel_scale=0.5,
    m_min=1e-3,
    m_max=1e-2,
):
    random.seed(seed)
    bodies = []
    for _ in range(N):
        bodies.append(
            Body(
                m=random.uniform(m_min, m_max),
                x=random.uniform(-pos_scale, pos_scale),
                y=random.uniform(-pos_scale, pos_scale),
                z=random.uniform(-pos_scale, pos_scale),
                vx=random.uniform(-vel_scale, vel_scale),
                vy=random.uniform(-vel_scale, vel_scale),
                vz=random.uniform(-vel_scale, vel_scale),
            )
        )
    return bodies


def time_simulation(bodies, solver, dt, steps, softening):
    cfg = SimulationConfig(dt=dt, timesteps=steps, softening=softening)
    sim = Simulation(
        bodies=bodies,
        cfg=cfg,
        integrator=LeapfrogIntegrator(),
        solver=solver,
    )

    t0 = time.perf_counter()
    sim.run()
    t1 = time.perf_counter()

    return t1 - t0, sim


def extract_accuracy_metrics(sim):
    return {
        "energy": max(abs(x) for x in sim.energy_drift),
        "angular": max(abs(x) for x in sim.angular_momentum_drift),
        "com": max(abs(x) for x in sim.com_drift),
    }


#sweeping theta for testing 

def theta_sweep(
    N,
    thetas,
    dt=2e-3,
    steps=2000,
    softening=1e-3,
    repeats=3,
):
    results = []

    print(f"\nPhase 5 θ sweep at N={N}")
    print("θ | runtime (s) | max ΔE | max ΔL | max COM")

    for theta in thetas:
        runtimes = []
        metrics = None

        for _ in range(repeats):
            bodies = make_random_bodies(N=N, seed=42)

            t, sim = time_simulation(
                bodies=[Body(*b.asTuple()) for b in bodies],
                solver=BarnesHutSolver(theta=theta),
                dt=dt,
                steps=steps,
                softening=softening,
            )

            runtimes.append(t)
            metrics = extract_accuracy_metrics(sim)

        result = {
            "theta": theta,
            "runtime": min(runtimes),
            "energy": metrics["energy"],
            "angular": metrics["angular"],
            "com": metrics["com"],
        }

        results.append(result)

        print(
            f"{theta:>3} | "
            f"{result['runtime']:.3f} | "
            f"{result['energy']:.2e} | "
            f"{result['angular']:.2e} | "
            f"{result['com']:.2e}"
        )

    return results


#plotting

def plot_theta_results(results):
    thetas = [r["theta"] for r in results]

    plt.figure()
    plt.plot(thetas, [r["runtime"] for r in results], marker="o")
    plt.xlabel("θ")
    plt.ylabel("Runtime (s)")
    plt.title("Barnes–Hut Runtime vs θ (N=1000)")
    plt.show()

    plt.figure()
    plt.plot(thetas, [r["energy"] for r in results], marker="o")
    plt.yscale("log")
    plt.xlabel("θ")
    plt.ylabel("Max Relative Energy Drift")
    plt.title("Energy Drift vs θ (N=1000)")
    plt.show()

    plt.figure()
    plt.plot(thetas, [r["angular"] for r in results], marker="o")
    plt.yscale("log")
    plt.xlabel("θ")
    plt.ylabel("Max Angular Momentum Drift")
    plt.title("Angular Momentum Drift vs θ (N=1000)")
    plt.show()

    plt.figure()
    plt.plot(thetas, [r["com"] for r in results], marker="o")
    plt.yscale("log")
    plt.xlabel("θ")
    plt.ylabel("Max COM Drift")
    plt.title("Center-of-Mass Drift vs θ (N=1000)")
    plt.show()


#only runs when executed directly
if __name__ == "__main__":
    thetas = [0.3, 0.5, 0.7, 1.0]
    results = theta_sweep(N=1000, thetas=thetas)
    plot_theta_results(results)