####################################################################
# Imports and given constants and classes to be used in your code
# DO NOT CHANGE ANY OF THE BELOW
####################################################################

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle 
import math

# UNITS USED
# Time: Astronomical Units (AU). 1AU = distance between Sun and earth
# Mass: Solar Mass (M_sun). 1M_sun = mass of the Sun
# Time: Year (yr). 1yr = 1 year (one rotation of the Earth around the Sun)
# Luminosity: Solar Luminosity (L_sun). 1L_sun = luminosity of the Sun
# Velocity: AU/yr

G = 4*(math.pi)**2  # Universal gravitational constant (AU^3/M_sun/yr^2)
scoeff = 22.62      # Scorching coefficient, approx. 0.745 AU from the Sun (L_sun/AU^2)
fcoeff = 3.14       # Freezing coefficient, approx. 2 AU from the Sun (L_sun/AU^2)

class Box:
    def __init__(self,x0,y0,x1,y1):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
        self.maxSide = max(self.x1-self.x0,self.y1-self.y0)
        self.mx = (self.x0+self.x1)/2
        self.my = (self.y0+self.y1)/2
        
    def isIn(self,p):
        return self.x0 <= p.x <= self.x1 and self.y0 <= p.y <= self.y1
    
    def asTuple(self):
        return (self.x0,self.y0,self.x1,self.y1,self.mx,self.my)
    
    def split4(self):
        bNE = Box(self.mx,self.my,self.x1,self.y1)
        bNW = Box(self.x0,self.my,self.mx,self.y1)
        bSW = Box(self.x0,self.y0,self.mx,self.my)
        bSE = Box(self.mx,self.y0,self.x1,self.my)
        return (bNE,bNW,bSW,bSE)

    # Build a box that encloses all Bodies from an array 
    @staticmethod
    def getBox(P):
        if len(P) == 0: return None
        x0 = x1 = P[0].x
        y0 = y1 = P[0].y
        for p in P: 
            if p.x < x0: x0 = p.x
            if p.y < y0: y0 = p.y
            if p.x > x1: x1 = p.x
            if p.y > y1: y1 = p.y
        return Box(x0,y0,x1,y1)    
    
    def __str__(self):
        return f"Box({self.x0},{self.y0},{self.x1},{self.y1})"
        
class GNode:
    def __init__(self,box):
        self.box = box
        self.COM = None           # this node’s COM
        self.nbodies = 0          # number of bodies contained in node
        self.p = None             # if this node is a leaf with a Body
        self.children = None      # children are: [NE, NW, SW, SE]
        self.updateCOM()

    def isLeaf(self):
        return self.nbodies < 2
        
    def updateCOM(self):
        if self.isLeaf(): 
            ### if its a leaf node with no bodies
            if self.p == None:
                self.COM = Body(0, self.box.mx, self.box.my)
            ### otherwise calculate COM from the single body
            else:
                self.COM = Body(self.p.m, self.p.x, self.p.y)
            return
        # otherwise calculate COM from the 4 children
        x = y = m = 0
        for c in self.children:
            x += c.COM.x*c.COM.m
            y += c.COM.y*c.COM.m
            m += c.COM.m
        self.COM = Body(m, x/m, y/m)

    def niceStr(self): 
        S = ("├","─","└","│")
        angle = S[2]+S[1]+" "
        vdash = S[0]+S[1]+" "
        
        def niceRec(ptr,acc,pre,A):
            if ptr == None: raise Exception("A None GNode was found")
            val = f"{len(A)}:{ptr.box},{ptr.nbodies}"
            A.append(f"({ptr.COM.m}, {ptr.COM.x}, {ptr.COM.y})")
            if ptr.children==None: return acc+pre+val
            if pre == vdash: pre2 = S[3]+"  "
            elif pre == angle: pre2 = "   "
            else: pre2 = ""
            T = [vdash,vdash,vdash,angle]
            for i in range(4):
                T[i] = niceRec(ptr.children[i],acc+pre2,T[i],A)
            return acc+pre+val+"\n"+T[0]+"\n"+T[1]+"\n"+T[2]+"\n"+T[3]
            
        A = []
        s = niceRec(self,"","",A)+"\n"
        for i in range(len(A)):
            s += f"\n{i}{' '*(3-len(str(i)))}-> {A[i]}"
        return s
    
class Stack:
    def __init__(self):
        self.inList = []
        
    def push(self,v):
        self.inList.append(v)

    def pop(self):
        if len(self.inList) == 0: raise Exception("Popping from an empty stack")
        return self.inList.pop()
    
    def isEmpty(self):
        return len(self.inList) == 0
    
    def size(self):
        return len(self.inList)

    def toArray(self):
        return self.inList
    
    def __str__(self):
        return str(self.inList)





####################################################################
# Question 1: Intro to the n-body problem
# To solve: closest distance, predict good eras in 3-body problem
####################################################################

class Body:
    def __init__(self,m,x,y,vx=0,vy=0):
        self.m = m     # mass
        self.x = x     # position
        self.y = y     # position
        self.vx = vx   # velocity
        self.vy = vy   # velocity
        
    def squareDist(self,other):
        return (self.x-other.x)**2+(self.y-other.y)**2

    def __str__(self):
        return f"Body({self.m},{self.x},{self.y},{self.vx},{self.vy})"    

    def __repr__(self):
        return str(self)    

    def asTuple(self):
        return (self.m,self.x,self.y,self.vx,self.vy)    
    
    # Gives the next position and velocity of the current Body to its position
    # and velocity after time dt las elapsed, taking into account the pull
    # forces from the bodies in the array Bodies.
    def next(self, Bodies, dt):
        ret = Body(self.m, self.x, self.y)
        ax = ay = 0       
        # for each p in Bodies we compute their force on this Body and add it to its acceleration 
        for p in Bodies:
            if p == self: continue # current Body does not affect itself          
            # euclidian distance between p and this Body
            sq_distance = self.squareDist(p)
            # see e.g. https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation
            # for the vector form of Newton's law of gravity            
            ax += (p.x - self.x) *  p.m * G / (sq_distance**1.5)
            ay += (p.y - self.y) *  p.m * G / (sq_distance**1.5)
        # we compute displacement due to acceleration that we computed above
        # and due to its current speed because of inertia        
        ret.x += dt*dt*ax + dt*self.vx
        ret.y += dt*dt*ay + dt*self.vy
        # compute the velocity vectors
        ret.vx = (ret.x-self.x)/dt
        ret.vy = (ret.y-self.y)/dt

        return ret

    # Predict years of stability that self planet is going to have in a 3-sun solar system.
    # It suffices to take into account just the three suns (sunA,sunB,sunC) and the planet.
    # Stability is broken if one of the following criteria is violated: 
    # a. No scorching: lA/dA + lB/dB + lC/dC must be less than scoeff
    # b. No freezing: lA/dA + lB/dB + lC/dC must be greater than fcoeff
    # Notes: 
    # - dA, dB, dC are the squared distances between sunA, sunB, sunC and self, and 
    # - lA, lB, lC are the luminosities of sunA, sunB, sunC 
    # respectively. 
    def threeBodyProblem(self,sunA,sunB,sunC,lA,lB,lC): # 15%
        Newsim = Simulation([sunA,sunB,sunC,self])
        habitable = 0
        pss = Newsim.run()

        scoeff = 22.62      # Scorching coefficient, approx. 0.745 AU from the Sun (L_sun/AU^2)
        fcoeff = 3.14       # Freezing coefficient, approx. 2 AU from the Sun (L_sun/AU^2)

        for t in range(Newsim.timesteps):
            illum = (
                lA / pss[t][3].squareDist(pss[t][0]) +
                lB / pss[t][3].squareDist(pss[t][1]) +
                lC / pss[t][3].squareDist(pss[t][2])
            )



        ## if planet isnt scorched or frozen, add time step to habitable time
            if fcoeff < illum < scoeff:
                habitable += Newsim.dt
            else:
                return t/100 ### to get years

        return Newsim.total_time + Newsim.dt ### if never uninhabitable, return total time simulated + timestep
    
class Simulation:
    # Default simulation time 10yr, step time 0.01yr
    def __init__(self, Bodies, total_time = 10, dt=0.01):
        self.bodies = Bodies
        self.total_time = total_time
        self.dt = dt
        self.timesteps = int(total_time/dt)
        
    # Runs the simulation and produces an array of arrays of Bodies.
    # The t-th entry in the list is the position of Bodies after the t-th timestep
    # has been simulated    
    def run(self):
        pss = [None]*(self.timesteps+1)
        pss[0] = self.bodies
        for t in range(self.timesteps):
            # for every Body in the current timestep add its next position in next timestep
            pss[t+1] = [pss[t][i].next(pss[t],self.dt) for i in range(len(self.bodies))]
        return pss
    
    def closestDistance(self): # 15%
        ret = None
        pss = self.run()
        findMin = lambda a, b: a if (a is not None) and (a < b) else b

        for step in pss:
            for body in step:   
                for other in step:
                    if body != other:
                        ret = findMin(ret, body.squareDist(other))

        return ret if ret is None else math.sqrt(ret)
     
    # Shows an animation of the simulation using matplotlib
    def show(self,x0,y0,x1,y1):
        pss = self.run()
        # get figure and axes objects
        fig, ax = plt.subplots()
        # set some reasonable zoom on axes
        ax.set_xlim(x0,x1)
        ax.set_ylim(y0,y1)
        # put labels on the window for Bodies with different colours
        scatter = []
        for i in range(len(self.bodies)):
            scatter.append(ax.scatter([], [], marker='o', label=f'Body {i}'))
         # add timestep text to the legend
        time_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, 
                       verticalalignment='top', bbox=dict(boxstyle='round', 
                       facecolor='wheat', alpha=0.5))
        # update function from one position to the next
        def update(frame):
            for i in range(len(self.bodies)):
                scatter[i].set_offsets([pss[frame][i].x, pss[frame][i].y])
            time_text.set_text(f'Timestep: {frame}')                
            return scatter + [time_text]
        # generate and show the animation
        a = FuncAnimation(fig, update, frames=len(pss), interval=1, blit=True, repeat=False)
        plt.xlabel("X coordinate (AU)")
        plt.ylabel("Y coordinate (AU)")
        plt.title("Celestial body trajectories")
        plt.legend()
        plt.show()


####################################################################
# Question 2: The tree Gadget
# To solve: add and remove elements on a Gadget
####################################################################

class Gadget:
    def __init__(self,box):
        self.root = GNode(box)
        self.size = 0
        
    # Add Body p to this Gadget [30% TODO]        
    def add(self,p): 
        self._addingNode(self.root, p)
        self.size += 1
        return



### if a child has changed, we MUST update the COM of all its parents (as well as the parent's parent and etc)
### While descending the tree to insert p, we need to keep a stack of all nodes we visit.
### After insertion,go back up the stack (from leaf to root) and call updateCOM() on each node
    def _addingNode(self, node, p):

        parents = Stack()
        while True:
            # increase count for this node
            node.nbodies += 1
            parents.push(node) # add it to the stack

            # If it's an internal node, descend further down the gadget
            if node.children is not None:
                node = self.helper_add(node, p)
                continue

            # if it's a leaf and empty, place p
            elif node.p is None:
                node.p = p
                node.updateCOM()
                break

            else:
                # if it's a leaf but already has a body, split
                p0 = node.p ### the body already in the node
                node.p = None #### node is no longer a leaf (as we need to split) so this attrbute is cleared

                # create children
                NE, NW, SW, SE = node.box.split4()
                node.children = [GNode(NE), GNode(NW), GNode(SW), GNode(SE)]

                # place the OLD body
                child0 = self.helper_add(node, p0)
                child0.p = p0 ### restore attrbute p to make it a leaf node with 1 node
                child0.nbodies = 1 ### update nbodies
                child0.updateCOM()

                # now place the NEW body
                node = self.helper_add(node, p)
                # loop continues; if node is leaf&empty it will add p immediately


        while not parents.isEmpty():
            n = parents.pop()
            n.updateCOM()



    def helper_add(self, node, p):
        mx = node.box.mx
        my = node.box.my

        # Deterministic quadrant selection so that midline bodies are always in NE, as in spec.
        # NE
        if p.x >= mx and p.y >= my:
            return node.children[0]
        # NW
        if p.x <= mx and p.y >= my:
            return node.children[1]
        # SW
        if p.x <= mx and p.y <= my:
            return node.children[2]
        # SE
        return node.children[3]


    
    # Remove Body p from this Gadget [20% TODO]

    #the following methods work together to remove a body from the Gadget while updating all COMs and nbodies of parents.
    def remove(self, p): 
        if self.size == 0:
            return

        node = self.root
        parents = Stack()

        # DESCEND to the leaf that should contain p
        while True:
            parents.push(node)

            if node.children is None:
                break

            mx = node.box.mx
            my = node.box.my

            # choose the correct child index deterministically
            if p.x >= mx and p.y >= my:      # NE
                idx = 0
            elif p.x <= mx and p.y >= my:    # NW
                idx = 1
            elif p.x <= mx and p.y <= my:    # SW
                idx = 2
            else:                            # SE
                idx = 3

            child = node.children[idx]

            # only descend if subtree has bodies
            if child.nbodies == 0:
                return  # p not present

            node = child

        # verify leaf actually contains p (manual comparison, no built-ins)
        if node.p is None:
            return

        ### essentially does the same thing as the node.p.asTuple() built in
        ### without using the built in method
        if not (
            node.p.m  == p.m  and
            node.p.x  == p.x  and
            node.p.y  == p.y  and
            node.p.vx == p.vx and
            node.p.vy == p.vy
        ):
            return

        # remove p
        node.p = None

        # update all ancestors
        self._updateParents(parents)
        self.size -= 1


    def _updateParents(self, parents):
        
        while not parents.isEmpty(): #will continue until all parents are visited
            node = parents.pop()   # start from leaf going upwards (LIFO)
            node.nbodies -= 1 

            #if node is empty the it becomes a leaf
            if node.nbodies == 0: 
                node.p = None
                node.children = None
                node.updateCOM()
                continue 
    
            # if node has exactly one body left and is internal, we collapse it
            if node.nbodies == 1 and node.children is not None:
                # collapse to leaf
                leaf = self._findLeaf(node)
                node.children = None
                node.p = leaf
                node.updateCOM()
                continue
    
            node.updateCOM() #update the COM before next Node



    #used for collapsing the children
    '''def _checkingChildren(self, node):
        if node.children is None: #returns the leaf (stop case)
            return node.p
        for c in node.children:
            if c.nbodies > 0:
                return self._checkingChildren(c) #recursively checks for the leaf'''
    

    def _findLeaf(self, node):
        if node.children is None:
            return node.p

        for c in node.children:
            if c.nbodies > 0:
                return self._findLeaf(c)



    

    
    def __str__(self):
        return self.root.niceStr()

    # Collect all Bodies from this Gadget and return them in an array
    def getBodies(self):
        A = [None]*self.size
        i = 0; stack = Stack(); stack.push(self.root)
        while not stack.isEmpty():
            n = stack.pop()
            if n.isLeaf():
                if n.p is not None:
                    A[i] = n.p
                    i += 1
            else:
                for c in n.children: stack.push(c)
        return A   

    # Build a new Gadget and add in it all Bodies from array ps
    @staticmethod
    def fromBodies(ps):
        # calculate bottom-left and top-right positions
        if ps == []: return None
        x0 = x1 = ps[0].x
        y0 = y1 = ps[0].y
        for p in ps:
            if p.x < x0: x0 = p.x
            elif p.x > x1: x1 = p.x
            if p.y < y0: y0 = p.y
            elif p.y > y1: y1 = p.y
        # build gadget and add bodies 
        g = Gadget(Box(x0,y0,x1,y1))
        for p in ps: g.add(p)
        return g    
    
    @staticmethod
    def _drawNode(ax, node, showCom):
        x0, y0, x1, y1, mx, my = node.box.asTuple()
        w = x1-x0; h = y1-y0
        # draw the node rectangle
        rect = Rectangle((x0, y0), w, h, fill=False, linewidth=0.8, alpha=0.5)
        ax.add_patch(rect)
        # draw COM for internal nodes
        if showCom and not node.isLeaf():
            ax.plot([node.COM.x], [node.COM.y], marker='+', markersize=6)
            ax.plot([mx, node.COM.x], [my, node.COM.y], linewidth=0.3)            
        # draw leaf contents
        if node.isLeaf(): 
            if node.p is not None:
                ax.plot([node.p.x], [node.p.y], marker='o', markersize=3)
        else:
            # recurse into children [NE, NW, SW, SE]
            for c in node.children:
                if c is not None:
                    Gadget._drawNode(ax, c, showCom)
                
    def plot(self, figsize=(6,6)): # figsize : (w,h) in inches
        showCom=True         # Draw line with '+' from center to COM of internal nodes
        margin_ratio=0.05    # Extra margin around the root box
        fig, ax = plt.subplots(figsize=figsize)
        # Plot node bounds
        Gadget._drawNode(ax, self.root, showCom)
        # collect Bodies from leaves
        A = self.getBodies()
        if A != []:
            xs, ys = zip(*[(A[i].x,A[i].y) for i in range(len(A))])
            ax.plot(xs, ys, linestyle='none', marker='o', markersize=3)
        # Set view limits from root box (with margin)
        x0, y0, x1, y1, mx, my = self.root.box.asTuple()
        dx, dy = x1-x0, y1-y0
        pad_x = margin_ratio * max(dx, 1e-12)
        pad_y = margin_ratio * max(dy, 1e-12)
        ax.set_xlim(x0-pad_x, x1+pad_x)
        ax.set_ylim(y0-pad_y, y1+pad_y)
        ax.set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.show()



####################################################################
# Question 3: Using the Gadget to speed up simulations
# To solve: get bodies from a Gadget wrt reference body 
####################################################################

class FastSimulation(Simulation):
    # Runs the simulation and produces an array of arrays of Bodies.
    # The t-th entry in the list is the position of Bodies after the t-th timestep
    # has been simulated.
    # Simulation uses a Gadget in order to approximate for each body the list 
    # of other bodies gravitationally affecting its trajectory.
    def run(self,test=None): 
        pss = [None]*(self.timesteps+1)
        pss[0] = self.bodies
        for t in range(self.timesteps):
            # calculate the current Gadget
            g = Gadget.fromBodies(pss[t])
            # for every Body in the current timestep add its next position in next timestep
            # but using the gadget g
            A = pss[t][:]
            for i in range(len(A)):
                new_ps = FastSimulation.getBodies(g,A[i])
                A[i] = A[i].next(new_ps,self.dt)
                if test is not None: test[i] = new_ps
            pss[t+1] = A
        return pss

    # Barnes-Hut criterion deciding whether node n should be opened when calculating
    # the gravitational forces exerted on body p.
    # By defualt, if p is in n then we open.
    @staticmethod
    def BarnesHut(n,p):
        theta = 0.7
        if n.box.isIn(p): return True
        sqDist = (n.COM.x-p.x)**2+(n.COM.y-p.y)**2
        return n.box.maxSide**2/sqDist >= theta**2
        
    # Get bodies from gadget g in an optimal way, using criterion shouldOpen and body p to 
    # decide whether nodes in the gadget should be "opened" (i.e. their subnodes examined)
    # or not (i.e. the whole subtree approximated by its COM).
    @staticmethod
    def getBodies(g,p,shouldOpen=BarnesHut): # 20%
        ret = Stack()
        FastSimulation.getBodiesAux(g.root, p, ret, shouldOpen)
        ret = [ret.pop() for _ in range(ret.size())] #Changes it back to the original array
        return ret
    
    @staticmethod
    def getBodiesAux(n, p, bodies, shouldOpen=BarnesHut):
        if n is None:
            return
        
        # Base case to determine whether or not its a leaf.
        if n.nbodies <= 1:
            if n.p is not None:
                bodies.push(n.p)
            return 
        
        # Base case to determine whether or not the node shouldOpen (should continue to traverse)
        if not shouldOpen(n, p):
            bodies.push(n.COM)
            return
        
        # Now traverses each of the children gathering the approximated/actual bodies that affect p
        if n.children is None: return
        FastSimulation.getBodiesAux(n.children[0], p, bodies, shouldOpen) #NE
        FastSimulation.getBodiesAux(n.children[1], p, bodies, shouldOpen) #NW
        FastSimulation.getBodiesAux(n.children[2], p, bodies, shouldOpen) #SW
        FastSimulation.getBodiesAux(n.children[3], p, bodies, shouldOpen) #SE
        return
         




