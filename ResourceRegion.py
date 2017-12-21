import math


class ResourceRegion:
    def __init__(self, c, m):
        self.rootCity = c
        self.myMap = m
        self.culture = self.rootCity.culture
        self.nodes = []
        self.addNode(self.rootCity.node)
        self.resources = [0, 0]
        self.updateReg()

    def expungeReg(self):
        self.myMap.resourceRegions.remove(self)

    def addNode(self, node):
        if node.resourceRegion is not None:
            node.resourceRegion.nodes.remove(node)
        node.resourceRegion = self
        self.nodes.append(node)

    def cityCount(self):
        c = 0
        q = 0
        for p in self.nodes:
            if p.city is not None:
                c += 1
                q += p.city.population
        self.totalCities = c
        self.totalPop = q

    def sumResources(self):
        self.resources = [0, 0]  # [Food resources, industrial resources]
        rawPlant = 0
        rawMetal = 0
        rawAnimal = 0
        for p in self.nodes:
            if p.landmass is not None:
                rawPlant += p.vegetation
                rawMetal += p.metallicity
                rawAnimal += p.herbivores + (p.carnivores / 3)
            else:
                rawAnimal += 0.1
                rawPlant += 0.1
        m = self.culture.value.mainValues
        if "metallurgists" in m:
            rawMetal = rawMetal * 1.25
        if "agriculture" in m:
            rawPlant *= 1.15
        if "warriors" in m:
            rawAnimal *= 1.2
        if "simplicity" in m:
            rawMetal *= 0.75
            rawAnimal *= 0.75
            rawPlant *= 0.75
        scale = self.myMap.resourceScale
        self.resources[0] = scale * (rawPlant + (rawAnimal * 0.7))
        self.resources[1] = scale * (rawMetal + (rawAnimal * 0.3) + (rawPlant * 0.4))
        if "craftsmen" in m:
            self.resources[1] *= 1.1
        if "builders" in m:
            self.resources[1] *= 1.1
        if "collectivists" in m:
            self.resources[0] *= 1.05

    def updateReg(self):
        for p in self.nodes:
            if p.city is not None:
                p.resourceDist = math.log2(p.city.population)
        for p in self.nodes:
            for k in p.neighbors:
                if k.resourceRegion != self and k.resourceRegion is not None:
                    # In this case, it's a different region OTHER than noRegion.
                    if k.resourceRegion.culture != self.culture:
                        # Interact with a different region of a different society...
                        # This is a placeholder. In the future there will probably be wars here.
                        k.resourceRegion = k.resourceRegion
                    elif k.resourceRegion.culture == self.culture:
                        # Interacting with a different region of the same society...
                        if k.resourceDist < p.resourceDist:
                            self.addNode(k)
                            k.resourceDist = (p.resourceDist - 1) / 2
                elif k.resourceRegion is None:
                    if k.resourceDist < p.resourceDist:
                        self.addNode(k)
                        k.resourceDist = (p.resourceDist - 1) / 2
                elif k.resourceRegion == self:
                    if k.resourceDist < p.resourceDist:
                        k.resourceDist = (p.resourceDist - 1) / 2
        self.cityCount()
        self.sumResources()
        if len(self.nodes) == 0:
            self.expungeReg()
            return -1
