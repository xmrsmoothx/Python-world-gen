from Tools import *
from Node import Node
from River import River


class Landmass:
    def __init__(self, rootNode, sLevel):
        self.sealevel = sLevel
        self.color = self.landmassColor()
        self.root = rootNode
        self.nodes = []
        self.addNode(rootNode)
        self.boundary = []
        self.rivers = []
        self.fill()
        self.size = len(self.nodes)
        self.centermass()
        self.culturalNames = {}
        self.landmassType = self.lType(self.size)
        self.centroid = None

    def centermass(self):
        xTotal = sum([p.x for p in self.nodes])
        yTotal = sum([p.y for p in self.nodes])
        xx = xTotal / self.size
        yy = yTotal / self.size
        self.centroid = Node(xx, yy)

    def lType(self, s):
        if s <= 3:
            t = "isle"
        elif s <= 16:
            t = "atoll"
        elif s <= 512:
            t = "island"
        elif s <= 2048:
            t = "land"
        else:
            t = "continent"
        return t

    def landmassColor(self):
        h = math.floor(random.random() * 255)
        s = 128 + math.floor(random.random() * 128)
        v = 128 + math.floor(random.random() * 128)
        return (h, s, v)

    def addNode(self, p):
        if p not in self.nodes:
            self.nodes.append(p)
            p.landmass = self

    def addBoundary(self, p):
        if p not in self.boundary:
            self.boundary.append(p)

    def fill(self):
        for p in self.nodes:
            for k in p.neighbors:
                if k.elevation < self.sealevel:
                    self.addBoundary(p)
                while k not in self.nodes and k.elevation >= self.sealevel:
                    self.addNode(k)

    def addRiver(self, length):
        root = self.boundary[int(math.floor(random.random() * (len(self.boundary) - 1)))]
        riverLen = random.randint(length / 4, length)
        newRiver = River(root, riverLen, self)
        self.rivers.append(newRiver)

    def cullStreams(self, minLen):
        for r in self.rivers:
            if len(r.nodes) < minLen:
                self.rivers.remove(r)
                r.removeRiver()
