class bodyWater:
    def __init__(self, rootNode, sLevel):
        self.sealevel = sLevel
        self.root = rootNode
        self.nodes = []
        self.addNode(rootNode)
        self.fill()
        self.culturalNames = {}

    def addNode(self, p):
        if p not in self.nodes:
            self.nodes.append(p)
            p.bodyWater = self

    def fill(self):
        for p in self.nodes:
            for k in p.neighbors:
                while k not in self.nodes and k.elevation < self.sealevel:
                    self.addNode(k)
