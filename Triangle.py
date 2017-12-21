from Tools import *


class Triangle:
    def __init__(self, aa, bb, cc):
        self.verts = [None, None, None]
        self.verts[0] = aa
        self.verts[1] = bb
        self.verts[2] = cc
        self.neighbors = []

    def sharesNeighbors(self, other):
        shares = 0
        if self.verts[0] in other.verts:
            shares += 1
        if self.verts[1] in other.verts:
            shares += 1
        if self.verts[2] in other.verts:
            shares += 1
        return shares

    def toString(self):
        self.selfString = "("
        self.selfString += self.verts[0].toString()
        self.selfString += ", "
        self.selfString += self.verts[1].toString()
        self.selfString += ", "
        self.selfString += self.verts[2].toString()
        self.selfString += ")"
        return self.selfString

    def printSelf(self):
        print(self.toString())

    def drawGraph(self, drawer):
        drawer.line([self.verts[0].coords(), self.verts[1].coords()], "black", 1)
        drawer.line([self.verts[1].coords(), self.verts[2].coords()], "black", 1)
        drawer.line([self.verts[2].coords(), self.verts[0].coords()], "black", 1)

    def drawElevation(self, drawer, sealevel=0):
        elevationList = [self.verts[f].elevation for f in range(len(self.verts))]
        elevation = sum(elevationList) / 3
        col = math.floor(elevation * 255)
        dCol = (0, 0, col)
        underwater = 0
        for f in self.verts:
            if f.elevation < sealevel:
                underwater = 1
        if underwater == 1:
            dCol = (142, 128, clamp(col, 64, 76))
        drawer.polygon([self.verts[0].coords(), self.verts[1].coords(), self.verts[2].coords()], fill=dCol,
                       outline=dCol)

    def drawWildlife(self, drawer, sealevel=0):
        carnivores = math.floor(128 * sum([self.verts[f].carnivores for f in range(len(self.verts))]) / 3)
        herbivores = math.floor(128 * sum([self.verts[f].herbivores for f in range(len(self.verts))]) / 3)
        dCol = (128 + (carnivores) - (herbivores), (carnivores + herbivores) * 2, 128)
        underwater = 0
        for f in self.verts:
            if f.elevation < sealevel:
                underwater = 1
        if underwater == 1:
            dCol = (142, 128, 70)
        drawer.polygon([self.verts[0].coords(), self.verts[1].coords(), self.verts[2].coords()], fill=dCol,
                       outline=dCol)

    def drawLandmass(self, drawer, sealevel=0):
        landmass = None
        for p in self.verts:
            if p.landmass is not None:
                landmass = p.landmass
        if landmass is not None:
            dCol = landmass.color
        else:
            dCol = "black"
        underwater = 0
        for f in self.verts:
            if f.elevation < sealevel:
                underwater = 1
        if underwater == 1:
            dCol = (142, 64, 64)
        drawer.polygon([self.verts[0].coords(), self.verts[1].coords(), self.verts[2].coords()], fill=dCol,
                       outline=dCol)

    def drawReal(self, drawer, sl):
        elevationList = [self.verts[f].elevation for f in range(len(self.verts))]
        elevation = sum(elevationList) / 3
        col = math.floor(elevation * 255)
        underwater = 0
        avgHue = math.floor(sum([self.verts[f].biomeColor[0] for f in range(len(self.verts))]) / 3)
        avgSat = math.floor(sum([self.verts[f].biomeColor[1] for f in range(len(self.verts))]) / 3)
        avgVal = math.floor(sum([self.verts[f].biomeColor[2] for f in range(len(self.verts))]) / 3)
        dCol = (avgHue, avgSat, avgVal)
        for f in self.verts:
            if f.elevation < sl:
                underwater = 1
        if underwater == 1:
            dCol = (142, 128, clamp(col, 64, 76))
        drawer.polygon([self.verts[0].coords(), self.verts[1].coords(), self.verts[2].coords()], fill=dCol,
                       outline=dCol)
