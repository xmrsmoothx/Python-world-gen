class Region:
    def __init__(self, rootNode):
        self.root = rootNode
        self.nodes = []
        self.culturalNames = {}
        self.addNode(self.root)
        self.biome = self.root.biome
        self.fill()
        if self.biome == "boreal" or self.biome == "tropical":
            self.biome += " forest"
        if self.biome == "frozen":
            self.biome += " tundra"
        if self.biome == "water":
            bodysize = len(self.nodes)
            if bodysize > 1024:
                self.biome = "ocean"
            elif bodysize > 512:
                self.biome = "sea"
            else:
                self.biome = "lake"

    def addNode(self, p):
        if p not in self.nodes:
            self.nodes.append(p)
            p.region = self

    def fill(self):
        for p in self.nodes:
            for k in p.neighbors:
                if k.biome == self.biome:
                    self.addNode(k)
