
## Contains class for body and a SystemState used to store a list of bodies 


from typing import List
import math

# Gravitational constant in AU^3 / (M_sun * yr^2)
G = 4 * (math.pi)**2


class Body:
    def __init__(self, m, x, y, vx=0, vy=0):
        self.m = m     # mass
        self.x = x     # position x
        self.y = y     # position y
        self.vx = vx   # velocity x
        self.vy = vy   # velocity y

    def squareDist(self, other):
        return (self.x - other.x)**2 + (self.y - other.y)**2

    def __str__(self):
        return f"Body({self.m},{self.x},{self.y},{self.vx},{self.vy})"

    def __repr__(self):
        return str(self)

    def asTuple(self):
        return (self.m, self.x, self.y, self.vx, self.vy)


class SystemState:

    def __init__(self, bodies: List[Body]):
        self.bodies = bodies

    def copy(self):
        # deep copy of system state
        return SystemState([
            Body(b.m, b.x, b.y, b.vx, b.vy) 
            for b in self.bodies
        ])