from Tools import *


class Node:
    def __init__(self, xx, yy):
        self.biomeColor = None
        self.rainfall = None
        self.temp = None
        self.x = xx
        self.y = yy
        self.neighbors = []
        self.landmass = None
        self.bodyWater = None
        self.river = None
        self.region = None
        self.city = None
        self.culture = None
        self.allegiance = 0
        self.resourceRegion = None
        self.resourceDist = 0
        self.slope = None

    def coords(self):
        tupleVert = (self.x, self.y)
        return tupleVert

    def dist(self, n):
        distX = abs(self.x - n.x)
        distY = abs(self.y - n.y)
        dist = math.sqrt((distX ** 2) + (distY ** 2))
        return dist

    def isLinked(self, nNode):
        if nNode in self.neighbors:
            return 1
        else:
            return 0

    def link(self, newNeighbor):
        self.neighbors.append(newNeighbor)
        newNeighbor.neighbors.append(self)

    def hasWaterNeighbor(self, sealevel):
        for n in self.neighbors:
            if n.elevation < sealevel:
                return 1
        return 0

    def nearestNeighbor(self, n):
        distance = 100000000000
        nbr = self
        for p in self.neighbors:
            newDist = p.dist(n)
            if newDist < distance:
                nbr = p
                distance = newDist
        return nbr

    def westMostNeighbor(self):
        n = self
        for p in self.neighbors:
            if p.x < n.x:
                n = p
        self.westMostNbr = n
        return n

    def minNeighbor(self):
        n = self
        ne = self.elevation
        for p in self.neighbors:
            pe = p.elevation
            if pe < ne:
                ne = pe
                n = p
        return n

    def maxNeighbor(self):
        n = self
        ne = self.elevation
        for p in self.neighbors:
            pe = p.elevation
            if pe > ne:
                ne = pe
                n = p
        return n

    def getSlope(self):
        minimum = self.minNeighbor()
        maximum = self.maxNeighbor()
        rise = maximum.elevation - minimum.elevation
        run = self.dist(minimum) + self.dist(maximum)
        if run != 0:
            m = rise / run
        else:
            m = 0.1
        self.slope = m
        return m

    def slopeDirection(self):
        self.westMostNeighbor()
        if self.westMostNbr.elevation > self.elevation:
            self.slopeDir = -1
        else:
            self.slopeDir = 1

    def realSlope(self):
        self.getSlope()
        self.slopeDirection()
        self.rSlope = self.slope * self.slopeDir
        return self.rSlope

    def smooth(self):
        nbrs = []
        nbrs.append(self.elevation)
        for i in self.neighbors:
            nbrs.append(i.elevation)
        self.elevation = sum(nbrs) / len(nbrs)

    def setVegetation(self):
        tempFitness = 1 - abs(0.3 - self.temp)
        elevationFitness = 1 - abs(0.45 - self.elevation)
        fertilityFitness = self.fertility + (self.metallicity * 0.25)
        rainFitness = 1 - abs(self.rainfall - 0.8)
        vegFitness = ((tempFitness + 1) * (elevationFitness + 0.5) * (fertilityFitness + 0.5) * (rainFitness + 2))
        self.vegetation = clamp(vegFitness / 16, 0, 1)

    def setBiome(self, sl):
        if self.elevation > 0.9:
            self.biome = "mountain"
        elif self.temp < 0.2:
            if self.rainfall < 0.08:
                self.biome = "frozen"
            elif self.rainfall < 0.25:
                self.biome = "tundra"
            else:
                self.biome = "shrubland"
        elif self.temp < 0.3:
            if self.rainfall < 0.03:
                self.biome = "tundra"
            elif self.rainfall < 0.05:
                self.biome = "shrubland"
            elif self.rainfall < 0.3:
                self.biome = "boreal"
            elif self.rainfall < 0.4:
                self.biome = "forest"
            else:
                self.biome = "tropical"
        elif self.temp < 0.5:
            if self.rainfall < 0.04:
                self.biome = "desert"
            elif self.rainfall < 0.06:
                self.biome = "savanna"
            elif self.rainfall < 0.09:
                self.biome = "shrubland"
            elif self.rainfall < 0.2:
                self.biome = "forest"
            else:
                self.biome = "tropical"
        else:
            if self.rainfall < 0.05:
                self.biome = "desert"
            else:
                self.biome = "tropical"
        if self.elevation < sl:
            self.biome = "water"

    def claim(self, n):
        n.culture = self.culture
        n.allegiance = self.allegiance + 1
        if self.culture.name not in n.region.culturalNames:
            n.region.culturalNames[self.culture.name] = self.culture.language.genName()
        if n.landmass != None:
            if self.culture.name not in n.landmass.culturalNames:
                n.landmass.culturalNames[self.culture.name] = self.culture.language.genName()
        if n.river != None:
            if self.culture.name not in n.river.culturalNames:
                n.river.culturalNames[self.culture.name] = self.culture.language.genName()

    def updateAllegiance(self):
        if self.culture is not None:
            chance = 1 / clamp(self.allegiance, 0.00001, 512)
            for n in self.neighbors:
                roll = random.random()
                if roll <= chance:
                    self.claim(n)

    # NOTE: you can override the __str__() inbuilt method, which is the equivalent of what you are doing here
    # Advantage of overriding __str__():
    # You can write print(someNodeVariable) instead of print(someNodeVariable.toString())
    # As python automatically calls the __str__() class method when you try to print an object instance
    # This is why, when you do print(someObject) without overriding __str__(),
    #   you get something like <__main__.someObject object at 0x[RAM address]>
    def toString(self):
        self.selfString = "("
        self.selfString += str(self.x)
        self.selfString += ", "
        self.selfString += str(self.y)
        self.selfString += ")"
        return self.selfString

    def printSelf(self):
        print(self.toString())

    def drawPoint(self, drawer, radius, color):
        drawCircle(drawer, self.x, self.y, radius, color)

    def drawElevation(self, drawer, pts=0):
        col = math.floor(self.elevation * 255)
        dCol = (0, 0, col)
        if self.elevation >= 1:
            dCol = "red"
        if self.landmass is not None:
            for n in self.neighbors:
                if n.landmass == self.landmass:
                    drawer.line([self.coords(), n.coords()], dCol, 3)
        if pts == 1:
            self.drawPoint(drawer, 1, dCol)

    def drawLandmass(self, drawer, pts=0, drawNeighbors=0):
        if self.landmass is None:
            dCol = "black"
        else:
            dCol = self.landmass.color
        if self.landmass is not None:
            for n in self.neighbors:
                if n.landmass == self.landmass:
                    drawer.line([self.coords(), n.coords()], dCol, 3)
        if self.landmass is not None and drawNeighbors == 1:
            for n in self.neighbors:
                drawer.line([self.coords(), n.coords()], dCol, 1)
        if pts == 1:
            self.drawPoint(drawer, 1, "black")

    def drawReal(self, drawer, sl):
        dCol = self.biomeColor
        if self.landmass is not None and 0:
            for n in self.neighbors:
                if n.landmass == self.landmass:
                    drawer.line([self.coords(), n.coords()], dCol, 3)
