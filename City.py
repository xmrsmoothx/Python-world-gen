from Tools import *
from ResourceRegion import ResourceRegion
from Culture import Culture


class City:
    def __init__(self, n, pop=50, cltr=None, m=None):
        self.myMap = m
        self.myMap.cities.append(self)
        self.node = n
        self.node.city = self
        self.population = pop
        if cltr is None:
            self.culture = Culture(self.node, m)
        else:
            self.culture = cltr
        self.cityType = self.cType(self.population)
        self.name = self.culture.language.genName()
        self.region = self.node.region
        self.node.culture = self.culture
        self.node.allegiance = 1 / self.population
        if self.culture.name not in self.node.region.culturalNames:
            self.node.region.culturalNames[self.culture.name] = self.culture.language.genName()
        if self.culture.name not in self.node.landmass.culturalNames:
            self.node.landmass.culturalNames[self.culture.name] = self.culture.language.genName()
        for q in self.node.neighbors:
            if q.region.biome == "ocean" or q.region.biome == "lake" or q.region.biome == "sea":
                if self.culture.name not in q.region.culturalNames:
                    q.region.culturalNames[self.culture.name] = self.culture.language.genName()
        if self.node.river is not None:
            if self.culture.name not in self.node.river.culturalNames:
                self.node.river.culturalNames[self.culture.name] = self.culture.language.genName()
        reg = ResourceRegion(self, self.myMap)
        self.rawResources = [0, 0]

    def updateDemog(self):
        rscShare = self.population / self.node.resourceRegion.totalPop
        rscMax = [0, 0]
        rscMax[0] = rscShare * self.node.resourceRegion.resources[0]
        rscMax[1] = rscShare * self.node.resourceRegion.resources[1]
        mpo = 0.1  # Maximum personal output (max resources production per person)
        m = self.culture.value.mainValues
        if "agriculture" in m:
            mpo *= 1.1
        if "simplicity" in m:
            mpo *= 0.8
        if "warriors" in m:
            mpo *= 1.05
        if "builders" in m:
            mpo *= 1.05
        if "metallurgists" in m:
            mpo *= 1.05
        if "craftsmen" in m:
            mpo *= 1.05
        self.foodProduction = min(self.population * mpo, rscMax[0])
        self.industrialProduction = min(self.population * mpo, rscMax[1])
        mpc = 0.080  # Maximum personal consumption (max food needed per person)
        mpc -= 0.002 * (math.log10(clamp(self.population, 1, 1000000)) + 1)
        # As the population grows, need less food per person due to economies of scale
        mpc += 0.001 * math.log2(clamp(len(self.node.resourceRegion.nodes) - 16, 1, 1000000))
        # As the resource region grows, need more food per person due to transporation distance
        if "collectivists" in m:
            mpc -= 0.001
        if "freedom" in m:
            mpc += 0.001
        if "simplicity" in m:
            mpc -= 0.001
        if "greed" in m:
            mpc += 0.001
        self.foodConsumption = mpc * self.population
        diff = self.foodProduction - self.foodConsumption
        growth = clamp(diff / mpc, -self.population * 0.1, self.population * 0.1)
        self.population = math.ceil(self.population * 0.99)  # Age related death
        self.population = clamp(math.floor(self.population + growth + random.choice([-1, 0, 0, 0, 0, 1])), 1, 1000000)
        self.cityType = self.cType(self.population)

    def cType(self, p):
        if p <= 35:
            c = random.choice(["bivouac", "camp", "camp", "encampment", "campsite"])
        elif p <= 200:
            c = random.choice(["village", "village", "hamlet"])
        elif p <= 1000:
            c = random.choice(["township", "settlement"])
        elif p <= 5000:
            c = "town"
        elif p <= 20000:
            c = "city"
        else:
            c = "metropolis"
        return c.capitalize()

    def cityInfo(self):
        n = self.name + " (" + self.cityType + ")\n"
        n += "Population: " + str(self.population)
        return n

    def cultureInfo(self):
        info = self.culture.information()
        return info

    def drawTent(self, drawer, xx, yy, col, out):
        p0 = (xx, yy - 6)
        p1 = (xx - 4, yy + 2)
        p2 = (xx + 4, yy + 2)
        drawer.polygon([p0, p1, p2], outline=out, fill=col)
        p0 = (xx, yy - 1)
        p1 = (xx - 1, yy + 2)
        p2 = (xx + 1, yy + 2)
        drawer.polygon([p0, p1, p2], outline=out, fill=out)
        p0 = (xx, yy - 6)
        p1 = (xx - 2, yy - 8)
        p2 = (xx + 2, yy - 8)
        drawer.line([p0, p1], fill=out, width=1)
        drawer.line([p0, p2], fill=out, width=1)

    def drawHut(self, drawer, xx, yy, col, out):
        p0 = (xx + 3, yy + 2)
        p1 = (xx - 3, yy + 2)
        p2 = (xx - 3, yy - 3)
        p3 = (xx, yy - 6)
        p4 = (xx + 3, yy - 3)
        drawer.polygon([p0, p1, p2, p3, p4], outline=out, fill=col)
        p1 = (xx - 1, yy - 1)
        p2 = (xx + 1, yy + 2)
        drawer.rectangle([p1, p2], outline=out, fill=out)
        p0 = (xx, yy - 6)
        p1 = (xx + 5, yy - 1)
        p2 = (xx - 5, yy - 1)
        drawer.line([p0, p1], fill=out, width=1)
        drawer.line([p0, p2], fill=out, width=1)

    def drawVillage(self, drawer, xx, yy, col, out):
        p0 = (xx, yy - 5)
        p1 = (xx - 2, yy - 7)
        p2 = (xx - 4, yy - 5)
        p3 = (xx - 4, yy + 3)
        p4 = (xx, yy + 3)
        drawer.polygon([p0, p1, p2, p3, p4], outline=out, fill=col)
        p0 = (xx - 6, yy - 3)
        p2 = (xx + 2, yy - 3)
        drawer.line([p1, p0], fill=out, width=1)
        drawer.line([p1, p2], fill=out, width=1)
        p3 = (xx - 2, yy - 4)
        p4 = (xx - 2, yy - 3)
        drawer.line([p4, p3], fill=out, width=1)
        p0 = (xx, yy - 1)
        p1 = (xx + 2, yy - 3)
        p2 = (xx + 4, yy - 1)
        p3 = (xx + 4, yy + 3)
        p4 = (xx, yy + 3)
        drawer.polygon([p0, p1, p2, p3, p4], outline=out, fill=col)
        p0 = (xx - 2, yy + 1)
        p2 = (xx + 6, yy + 1)
        drawer.line([p1, p0], fill=out, width=1)
        drawer.line([p1, p2], fill=out, width=1)
        p3 = (xx + 2, yy + 3)
        p4 = (xx + 2, yy + 1)
        drawer.line([p4, p3], fill=out, width=1)

    def drawTown(self, drawer, xx, yy, col, out):
        p0 = (xx - 4, yy + 3)
        p1 = (xx - 4, yy - 2)
        p2 = (xx, yy - 5)
        p3 = (xx + 4, yy - 2)
        p4 = (xx + 4, yy + 3)
        drawer.polygon([p0, p1, p3, p4], outline=out, fill=col)
        drawer.polygon([p1, p2, p3], outline=out, fill=out)
        p0 = (xx - 2, yy + 4)
        p1 = (xx - 2, yy - 5)
        p2 = (xx, yy - 7)
        p3 = (xx + 2, yy - 5)
        p4 = (xx + 2, yy + 4)
        drawer.polygon([p0, p1, p3, p4], outline=out, fill=col)
        drawer.polygon([p1, p2, p3], outline=out, fill=out)
        p0 = (xx, yy)
        p1 = (xx, yy + 3)
        drawer.line([p1, p0], fill=out, width=1)

    def drawCity(self, drawer, xx, yy, col, out):
        p0 = (xx - 5, yy - 3)
        p1 = (xx - 5, yy + 1)
        p2 = (xx - 3, yy + 1)
        p3 = (xx - 3, yy - 5)
        drawer.polygon([p0, p1, p2, p3], outline=out, fill=col)
        p0 = (xx + 5, yy - 4)
        p1 = (xx + 5, yy + 1)
        p2 = (xx + 3, yy + 1)
        p3 = (xx + 3, yy - 4)
        drawer.polygon([p0, p1, p2, p3], outline=out, fill=col)
        p0 = (xx - 1, yy - 7)
        p1 = (xx - 1, yy + 2)
        p2 = (xx + 1, yy + 2)
        p3 = (xx + 1, yy - 7)
        drawer.polygon([p0, p1, p2, p3], outline=out, fill=col)
        p0 = (xx - 6, yy + 2)
        p1 = (xx - 6, yy - 1)
        p2 = (xx + 6, yy - 1)
        p3 = (xx + 6, yy + 2)
        drawer.polygon([p0, p1, p2, p3], outline=out, fill=col)
        p0 = (xx, yy)
        p1 = (xx, yy + 1)
        drawer.line([p1, p0], fill=out, width=1)

    def drawSelf(self, drawer):
        col = self.culture.bannerColor
        out = (0, 0, 0)
        if self.population <= 35:
            self.drawTent(drawer, self.node.x, self.node.y, col, out)
        elif self.population <= 200:
            self.drawHut(drawer, self.node.x, self.node.y, col, out)
        elif self.population <= 1000:
            self.drawVillage(drawer, self.node.x, self.node.y, col, out)
        elif self.population <= 5000:
            self.drawTown(drawer, self.node.x, self.node.y, col, out)
        elif self.population <= 20000:
            self.drawCity(drawer, self.node.x, self.node.y, col, out)
