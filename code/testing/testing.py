from code.nbody.bodies import Body
from code.nbody.engine import Simulation, SimulationConfig
from code.nbody.integrators.leapfrog import LeapfrogIntegrator
import matplotlib.pyplot as plt
import math

# -------------------------------
# Helper: Compute COM-frame Angular Momentum
# -------------------------------
def compute_com_angular_momentum(bodies):
    # total mass
    M = sum(b.m for b in bodies)

    # center of mass
    x_cm = sum(b.m * b.x for b in bodies) / M
    y_cm = sum(b.m * b.y for b in bodies) / M

    # center-of-mass velocity
    vx_cm = sum(b.m * b.vx for b in bodies) / M
    vy_cm = sum(b.m * b.vy for b in bodies) / M

    # compute L in COM frame
    L = 0.0
    for b in bodies:
        x = b.x - x_cm
        y = b.y - y_cm
        vx = b.vx - vx_cm
        vy = b.vy - vy_cm
        L += b.m * (x * vy - y * vx)

    return L


# -------------------------------
# Initial Conditions
# -------------------------------
sun = Body(1.0, 0, 0, 0, -0.001 * 2 * math.pi)
planet = Body(0.001, 1, 0, 0, 2 * math.pi)

cfg = SimulationConfig(dt=0.0002, timesteps=5000, softening=0.0)

bodies = [sun, planet]

sim = Simulation(bodies, cfg, integrator=LeapfrogIntegrator())
_ = sim.run()

# -------------------------------
# Compute COM-Frame Angular Momentum Drift
# -------------------------------
L_values = []
for system_state in sim.state_history:
    L_values.append(compute_com_angular_momentum(system_state.bodies))

L0 = L_values[0]
L_drift = [(L - L0) / abs(L0) for L in L_values]


# -------------------------------
# PLOTS
# -------------------------------

# 1. Total Energy
plt.figure()
plt.plot(sim.energy_history)
plt.title("Total Energy vs Time")
plt.xlabel("Timestep")
plt.ylabel("Energy")


# 2. Angular Momentum (COM-frame)
plt.figure()
plt.plot(L_values)
plt.title("Angular Momentum (COM Frame)")
plt.xlabel("Timestep")
plt.ylabel("L")


# 3. Angular Momentum Drift
plt.figure()
plt.plot(L_drift)
plt.title("Angular Momentum Drift (COM Frame)")
plt.xlabel("Timestep")
plt.ylabel("(L - L0) / |L0|" )
plt.show()

plt.figure()
plt.plot(sim.energy_drift)
plt.title("Relative Energy Drift")
plt.xlabel("Timestep")
plt.ylabel("(E - E0)/|E0|")

plt.figure()
plt.plot(sim.angular_momentum_drift)
plt.title("Relative Angular Momentum Drift (about origin)")
plt.xlabel("Timestep")
plt.ylabel("(L - L0)/|L0|")




