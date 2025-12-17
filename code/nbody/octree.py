from code.nbody.bodies import Body

class OctreeNode:
    def __init__(self, center, half_size):
        self.total_mass = 0.0
        self.center_of_mass = (0.0, 0.0, 0.0)
        self.half_size = half_size
        self.center = center    
        self.body = None
        self.children = None 


    def insert(self, body: Body):
        pass