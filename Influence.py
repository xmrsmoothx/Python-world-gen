class Influence:
    def __init__(self, mMap, myNode, root):
        self.node = myNode
        self.rootNode = root
        self.myMap = mMap
        self.influenceOutput = {}
        for o in self.myMap.influenceOutputs.keys():
            self.influenceOutput[o] = 0
        self.translateInfluence()

    def setOutput(self, influenceType, strength):
        for o in self.myMap.influenceOutputs:
            if o in self.myMap.influences[influenceType]:
                output = strength * self.myMap.influences[influenceType][o]
                self.influenceOutput[o] += output

    def translateInfluence(self):
        if self.rootNode == 1:
            for p in self.node.neighbors:
                ni = Influence(self.myMap, p, 0)
                for j in ni.influenceOutput.keys():
                    self.influenceOutput[j] += ni.influenceOutput[j] / len(self.node.neighbors)
        strength = 1
        inf = self.node.biome
        self.setOutput(inf, strength)
        strength = self.node.herbivores
        inf = "carnivores"
        self.setOutput(inf, strength)
        strength = self.node.carnivores
        inf = "herbivores"
        self.setOutput(inf, strength)
        strength = self.node.temp
        inf = "temperature"
        self.setOutput(inf, strength)
        strength = self.node.elevation
        inf = "elevation"
        self.setOutput(inf, strength)
        strength = self.node.slope
        inf = "slope"
        self.setOutput(inf, strength)
        strength = self.node.vegetation
        inf = "vegetation"
        self.setOutput(inf, strength)
