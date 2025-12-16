
## This computes the gravitational accelerations + energies on each body due to all others (further physics added here)


import math
from typing import List
from code.nbody.bodies import Body, G


def compute_accelerations(bodies: List[Body], cfg):
    N = len(bodies)

    ax = [0.0] * N
    ay = [0.0] * N
    az = [0.0] * N

    soft2 = cfg.softening * cfg.softening

    for i in range(N):
        bi = bodies[i]
        for j in range(i + 1, N):
            bj = bodies[j]

            dx = bj.x - bi.x
            dy = bj.y - bi.y
            dz = bj.z - bi.z
            r2 = dx * dx + dy * dy + dz * dz + soft2
            r = r2 ** 0.5
            r3 = r2 * r

            # force magnitude per unit mass
            f = G / r3

            ax_i = f * bj.m * dx
            ay_i = f * bj.m * dy
            az_i = f * bj.m * dz

            ax_j = -f * bi.m * dx
            ay_j = -f * bi.m * dy
            az_j = -f * bi.m * dz

            ax[i] += ax_i
            ay[i] += ay_i
            az[i] += az_i

            ax[j] += ax_j
            ay[j] += ay_j
            az[j] += az_j

    return ax, ay, az



def compute_kinetic_energy(bodies: List[Body]):
    Total = 0.0
    for b in bodies:
        Total += 0.5 * b.m * (b.vx ** 2 + b.vy ** 2 + b.vz **2)
    return Total

def compute_potential_energy(bodies: List[Body],cfg):
    total = 0.0
    length = len(bodies)

    for i in range(length):
        bi = bodies[i]
        for j in range(i + 1, length):
            bj = bodies[j]
            dx = bj.x - bi.x
            dy = bj.y - bi.y
            dz = bj.z - bi.z

            dist = math.sqrt(dx*dx + dy*dy + dz*dz + cfg.softening * cfg.softening)

            total += -G * bi.m * bj.m / dist
    return total


def compute_angular_momentum(bodies):
    total = 0.0
    for b in bodies:
        total += b.m * (b.x * b.vy - b.y * b.vx)
    return total

def compute_linear_momentum(bodies):
    px_total = 0.0
    py_total = 0.0
    for b in bodies:
        px_total += b.m * b.vx
        py_total += b.m * b.vy
    return (px_total, py_total)


def compute_center_of_mass(bodies):
    total_mass = 0.0
    x_cm = 0.0
    y_cm = 0.0
    for b in bodies:
        total_mass += b.m
        x_cm += b.m * b.x
        y_cm += b.m * b.y
    if total_mass > 0:
        x_cm /= total_mass
        y_cm /= total_mass
    return (x_cm, y_cm)



