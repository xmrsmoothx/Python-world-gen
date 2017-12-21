from Tools import *
from tkinter import *
from Node import Node
from noiseMaker import noiseMaker
from Landmass import Landmass
from bodyWater import bodyWater
from Region import Region
from City import City
from GenConfiguration import *
from PIL import Image, ImageDraw, ImageTk


class Map:
    def __init__(self, aAtlas, numNodes, mapDimX, mapDimY):
        self.triangles = None
        self.atlas = aAtlas
        self.n = numNodes
        self.xDim = mapDimX
        self.yDim = mapDimY
        self.landmasses = []
        self.waterBodies = []
        self.regions = []
        self.cities = []
        self.cultures = []
        self.resourceRegions = []
        self.resourceScale = 100
        self.sealevel = 0.4
        self.setNorth()
        self.biomeColors()
        self.displayNo = None
        self.infoGui = None

    def setNorth(self):
        nx = math.floor(random.random() * self.xDim)
        ny = math.floor(random.random() * self.yDim)
        side = random.randint(0, 3)
        if side == 0:
            nx = 0
        elif side == 1:
            nx = self.xDim
        elif side == 2:
            ny = 0
        else:
            ny = self.yDim
        self.north = Node(nx, ny)

    def nodeLat(self, n):
        return "Latitude: " + str(round(n.y / self.distScale))

    def nodeLong(self, n):
        return "Longitude: " + str(round(n.x / self.distScale))

    def nodeElevation(self, n):
        return "Elevation: " + str(round((n.elevation - self.sealevel) * self.eScale, 1)) + "m"

    def nodeTemp(self, n):
        if n.biome != "water":
            return "Temperature: " + str(round((n.temp * self.tempScale) - 30, 1)) + " degrees"
        else:
            return "Temperature: " + str(round((n.temp * self.tempScale * 0.3), 1)) + " degrees"

    def nodeRain(self, n):
        if n.biome != "water":
            rain = "Rainfall: " + str(round(((n.rainfall) ** 2) * self.rainfallScale, 1)) + "cm/yr"
        else:
            rain = "Rainfall: " + str(round((n.rainfall) * self.rainfallScale * 0.03, 1)) + "cm/yr"
        return rain

    def nodeBiome(self, n):
        return n.biome

    def nodeRegion(self, n):
        if n.region.culturalNames == {}:
            return "Unnamed " + n.region.biome + "\n"
        else:
            names = ""
            for f in n.region.culturalNames.keys():
                q = n.region.culturalNames[f] + " "
                q += n.region.biome
                q += " (" + f + " language)"
                names += q + "\n"
            return names

    def nodeLandmass(self, n):
        if n.landmass.culturalNames == {}:
            return "Unnamed " + n.landmass.landmassType + "\n"
        else:
            names = ""
            for f in n.landmass.culturalNames.keys():
                q = n.landmass.culturalNames[f] + " "
                q += n.landmass.landmassType
                q += " (" + f + " language)"
                names += q + "\n"
            return names

    def nodeRiver(self, n):
        if n.river.culturalNames == {}:
            return "Unnamed " + "river" + "\n"
        else:
            names = ""
            for f in n.river.culturalNames.keys():
                q = n.river.culturalNames[f] + " "
                q += "river"
                q += " (" + f + " language)"
                names += q + "\n"
            return names

    def nodeFertility(self, n):
        return "Soil fertility: " + str(math.floor(n.fertility * self.fertScale)) + "%"

    def nodeMetallicity(self, n):
        return "Ground metallicity: " + str(math.floor(n.metallicity * self.metalScale)) + "ppm"

    def nodeVegetation(self, n):
        return "Vegetation: " + str(math.floor(n.vegetation * self.vegScale)) + "p/km" + chr(0x00B2)

    def nodeHerbivores(self, n):
        return "Herbivores: " + str(math.floor(n.herbivores * self.wildlifeScale)) + "/km" + chr(0x00B2)

    def nodeCarnivores(self, n):
        return "Carnivores: " + str(math.floor(n.carnivores * self.wildlifeScale)) + "/km" + chr(0x00B2)

    def nodeCityInfo(self, n):
        cityInfo = n.city.cityInfo() + "\n" + "\n"
        cityInfo += n.city.cultureInfo()
        return cityInfo

    def nodeTerritory(self, n):
        territory = ""
        if n.culture is not None:
            territory = n.culture.shortName() + " territory"
        return territory

    def nodeResReg(self, n):
        reg = ""
        if n.resourceRegion is not None:
            reg = "Inside " + n.resourceRegion.rootCity.name + " outskirts" + "\n"
            reg += "Total food available: " + str(
                math.floor(n.resourceRegion.resources[0] * self.resourceScale)) + " t/year \n"
            reg += "Total industrial resources available: " + str(
                math.floor(n.resourceRegion.resources[1] * self.resourceScale)) + " t/year \n"
        return reg

    def infoScales(self):
        self.distScale = 12
        self.eScale = 2000
        self.tempScale = 105
        self.rainfallScale = 7628
        self.fertScale = 100
        self.metalScale = 140000
        self.vegScale = 10295
        self.wildlifeScale = 500
        self.resourceScale = 1

    def nodeInfo(self, n):
        self.divWidth = 48
        info = ""
        if n.city is not None:
            info += strDivider(self.divWidth) + "\n"
            info += self.nodeCityInfo(n) + "\n"
        if n.resourceRegion is not None:
            info += strDivider(self.divWidth) + "\n"
            info += self.nodeResReg(n) + "\n"
        info += self.nodeTerritory(n) + "\n"
        info += strDivider(self.divWidth) + "\n"
        info += self.nodeRegion(n) + "\n"
        if n.landmass is not None:
            info += self.nodeLandmass(n) + "\n"
        if n.river is not None:
            info += self.nodeRiver(n) + "\n"
        info += self.nodeElevation(n) + "\n"
        info += self.nodeTemp(n) + "\n"
        info += self.nodeRain(n) + "\n"
        if n.biome != "water":
            info += strDivider(self.divWidth) + "\n"
            info += self.nodeFertility(n) + "\n"
            info += self.nodeMetallicity(n) + "\n"
            info += self.nodeVegetation(n) + "\n"
            info += self.nodeHerbivores(n) + "\n"
            info += self.nodeCarnivores(n) + "\n"
        return info

    def nearestNode(self, xx, yy):
        n = self.atlas[0]
        minDist = 1000000
        search = Node(xx, yy)
        for p in self.atlas:
            dist = p.dist(search)
            if dist < minDist:
                minDist = dist
                n = p
        return n

    def nearestCity(self, xx, yy):
        if len(self.cities) == 0:
            return self.atlas[0]
        n = self.cities[0]
        minDist = 1000000
        search = Node(xx, yy)
        for p in self.cities:
            dist = p.node.dist(search)
            if dist < minDist:
                minDist = dist
                n = p
        return n

    def setSeaLevel(self, n):
        self.sealevel = n

    def perlinElevation(self, octaves, scale=1):
        noiseThing = noiseMaker(math.floor(self.xDim / scale), math.floor(self.yDim / scale))
        for p in self.atlas:
            p.elevation = noiseThing.turbulence(math.floor(p.x), math.floor(p.y), 2 ** octaves)

    def flatten(self):
        for p in self.atlas:
            p.elevation = 0.5

    def elevationAdd(self, amount):
        for p in self.atlas:
            p.elevation = clamp(p.elevation + amount, 0, 1)

    def randomizeElevation(self, anchor, rng):
        for p in self.atlas:
            p.elevation = anchor - (random.random() * (rng / 2))

    def clampElevation(self):
        for p in self.atlas:
            p.elevation = clamp(p.elevation, 0, 1)

    def smooth(self, strength):
        for l in range(strength):
            for p in self.atlas:
                p.smooth()

    def cullDots(self):
        for p in self.atlas:
            count = len(p.neighbors)
            for n in p.neighbors:
                if n.elevation < self.sealevel:
                    count -= 1
            if count <= 1:
                p = clamp(p.elevation, self.sealevel - 0.01, 255)

    def buildLandmass(self, root):
        if root.landmass is not None or root.elevation < self.sealevel:
            return -1
        root.landmass = Landmass(root, self.sealevel)
        self.landmasses.append(root.landmass)

    def buildAllLand(self):
        print("Building landmasses...")
        for p in self.atlas:
            self.buildLandmass(p)

    def buildBodyWater(self, root):
        if root.landmass is not None or root.bodyWater is not None:
            return -1
        root.bodyWater = bodyWater(root, self.sealevel)
        self.waterBodies.append(root.bodyWater)

    def buildAllWater(self):
        print("Building lakes/oceans...")
        for p in self.atlas:
            self.buildBodyWater(p)

    def buildRegions(self):
        print("Building world regions...")
        for p in self.atlas:
            if p.region is None:
                newReg = Region(p)
                self.regions.append(newReg)

    def addRiver(self, length):
        c = 0
        landmass = None
        while c < len(self.landmasses):
            landmass = self.landmasses[int(math.floor(random.random() * len(self.landmasses)))]
            if landmass.size > length * 8:
                c = 100000000
            c += 1
        landmass.addRiver(length)

    def addMinorRiver(self, count):
        for i in range(count):
            self.addRiver(self.n / 512)

    def addMajorRiver(self, count):
        for i in range(count):
            self.addRiver(self.n / 128)

    def cullStreams(self):
        for l in self.landmasses:
            l.cullStreams(self.n / 32)

    def addSineHill(self, xx, yy, maximum=0.25, radius=128):
        hillCenter = Node(xx, yy)
        for p in self.atlas:
            dist = p.dist(hillCenter)
            if dist <= radius:
                multiplier = distMod(dist, radius)
                if p.elevation > 0.75:
                    multiplier = multiplier * ((1 - p.elevation) ** 2)
                p.elevation = p.elevation + (maximum * multiplier)

    def addHill(self, xx, yy, maximum=0.25, radius=128):
        hillCenter = Node(xx, yy)
        for p in self.atlas:
            dist = p.dist(hillCenter)
            if dist <= radius:
                multiplier = random.uniform(0.99, 1.01)
                p.elevation = maximum * multiplier

    def addMountains(self, num=5, height=0.2):
        avgRad = self.xDim / 3.5
        for i in range(num):
            xx = math.floor(random.random() * self.xDim)
            yy = math.floor(random.random() * self.yDim)
            hillRad = avgRad * random.uniform(0.8, 1.25)
            hillHeight = height * random.uniform(0.8, 1.25)
            self.addSineHill(xx, yy, hillHeight, hillRad)

    def addHills(self, num=8, height=0.1):
        avgRad = self.xDim / 6
        for i in range(num):
            xx = math.floor(random.random() * self.xDim)
            yy = math.floor(random.random() * self.yDim)
            hillRad = avgRad * random.uniform(0.8, 1.25)
            hillHeight = height * random.uniform(0.8, 1.25)
            self.addSineHill(xx, yy, hillHeight, hillRad)

    def addShape(self, shape):
        if shape == "island" or shape == "volcanic":
            self.addSineHill(self.xDim / 2, self.yDim / 2, 0.4, radius=self.xDim / 1.5)
            if shape == "volcanic":
                self.smooth(4)
                self.elevationAdd(-0.1)
                self.addSineHill(self.xDim / 2, self.yDim / 2, 0.25, radius=self.xDim * 1.5)
                self.smooth(4)
                self.addHill(self.xDim / 2, self.yDim / 2, 0.45, radius=self.xDim / 11)
                self.addSineHill(self.xDim / 2, self.yDim / 2, 0. - 0.1, radius=self.xDim / 11)
        if shape == "shore" or shape == "highlands":
            corner = random.randint(0, 3)
            if corner == 0:
                xx = 0
                yy = 0
            elif corner == 1:
                xx = self.xDim
                yy = 0
            elif corner == 2:
                xx = self.xDim
                yy = self.yDim
            else:
                xx = 0
                yy = self.yDim
            self.addSineHill(xx, yy, 0.45, radius=self.xDim * 1.5)
            if shape == "highlands":
                self.addSineHill(xx, yy, 0.3, radius=self.xDim * 2)
                self.addMountains()
        if shape == "archipelago":
            self.addHills(16, 0.25)
        if shape == "plain":
            self.addSineHill(self.xDim / 2, self.yDim / 2, 0.4, radius=self.xDim * 5)
        if shape != "volcanic":
            self.addMountains()
            self.addHills()
            self.erode(3)
            self.smooth(3)

    def addRandomShape(self):
        shp = random.choice(["highlands", "plain", "volcanic", "shore", "archipelago", "island"])
        self.addShape(shp)

    def erode(self, strength=3):
        for j in range(strength):
            for p in self.atlas:
                if p.elevation < self.sealevel * 1.02 and p.elevation > self.sealevel * 0.9:
                    if p.hasWaterNeighbor(self.sealevel):
                        p.elevation = p.elevation * 0.94
                        for n in p.neighbors:
                            n.elevation = n.elevation * 0.98

    def nodeSlopes(self):
        for p in self.atlas:
            p.getSlope()

    def rainfall(self):
        for p in self.atlas:
            rainfall = 0.4 * random.uniform(0.95, 1.05)
            rainfall = ((rainfall * (0.4 - (p.temp + p.elevation))) + rainfall) / 2
            q = p
            distance = self.n / 256
            diff = 0
            for i in range(int(math.floor(distance))):
                q = p.westMostNeighbor()
                diff += q.elevation - p.elevation
            diff = diff / distance
            if diff > 0.2 and 0:
                rainfall = clamp(rainfall * (0.5 - diff), 0, 1)
            if p.river is not None:
                rainfall = rainfall * 1.5
            for l in p.neighbors:
                if l.elevation < self.sealevel:
                    rainfall *= 1.25
            p.rainfall = clamp(rainfall, 0, 1)

    def temperature(self):
        for p in self.atlas:
            latTemp = 0.375 * (p.dist(self.north) / self.xDim)
            elevationTemp = 0.575 * (1 - p.elevation)
            p.temp = clamp(elevationTemp + latTemp, 0, 1)

    def soilProperties(self):
        for p in self.atlas:
            p.metallicity = 0.5 * p.elevation * random.uniform(0.666, 1.5)
            if p.river is not None:
                p.metallicity *= 1.25
            p.metallicity = clamp(p.metallicity, 0, 1)
            fertilityBase = abs(p.elevation - (self.sealevel * 1.2)) * random.uniform(0.8, 1.25)
            if fertilityBase == 0:
                p.fertility = 1
            else:
                p.fertility = 0.07 / fertilityBase
            if p.river is not None:
                p.fertility *= 2
            p.fertility = clamp(p.fertility, 0, 1)

    def vegetation(self):
        self.temperature()
        self.soilProperties()
        self.rainfall()
        for p in self.atlas:
            p.setVegetation()

    def setBiomes(self):
        self.nodeSlopes()
        self.vegetation()
        for p in self.atlas:
            p.setBiome(self.sealevel)
            p.biomeColor = self.biomeColors[p.biome]
            slope = clamp(p.realSlope() * (2 ** 16), -48, 48)
            shade = math.floor((-16) + p.biomeColor[2] + slope + (((p.elevation + 1) ** 2) * 16))
            p.biomeColor = (p.biomeColor[0], p.biomeColor[1], shade)

    def setWildlife(self):
        for p in self.atlas:
            p.herbivores = clamp((p.vegetation * random.uniform(0.8, 1.25)) ** 2, 0, 1)
            p.carnivores = clamp(((p.herbivores * 1.5) ** 2), 0, 1)

    def biomeColors(self):
        bColors = {}
        bColors["desert"] = (16, 64, 142)
        bColors["savanna"] = (48, 64, 136)
        bColors["tundra"] = (128, 48, 142)
        bColors["shrubland"] = (72, 96, 128)
        bColors["boreal"] = (110, 86, 108)
        bColors["forest"] = (96, 128, 128)
        bColors["tropical"] = (78, 162, 96)
        bColors["frozen"] = (134, 32, 206)
        bColors["mountain"] = (134, 0, 136)
        bColors["water"] = (142, 64, 64)
        self.biomeColors = bColors

    def values(self):
        self.valuesOutputs = {"travelers": 0,
                              "craftsmen": 0,
                              "traders": 0,
                              "superstitious": 0,
                              "metallurgists": 0,
                              "worshippers": 0,
                              "freedom": 0,
                              "shamans": 0,
                              "astrologists": 0,
                              "materialists": 0,
                              "agriculture": 0,
                              "collectivists": 0,
                              "builders": 0,
                              "naturalists": 0,
                              "simplicity": 0,
                              "greed": 0,
                              "sailors": 0,
                              "warriors": 0}
        self.values = {}
        self.values["swimming"] = {"simplicity": 0.2,
                                   "sailors": 0.8,
                                   "astrologists": 0.15,
                                   "freedom": 0.35,
                                   "warriors": 0.05}
        self.values["food"] = {"agriculture": -0.1,
                               "greed": -0.25,
                               "materialists": 0.55,
                               "collectivists": 0.25,
                               "simplicity": 0.45,
                               "worshippers": -0.2,
                               "superstitious": -0.2,
                               "traders": 0.2,
                               "warriors": 0.2,
                               "builders": 0.25}
        self.values["darkness"] = {"travelers": 0.2,
                                   "collectivists": -0.4,
                                   "superstitious": 0.7,
                                   "greed": 0.4,
                                   "astrologists": 0.25,
                                   "materialists": -0.15,
                                   "shamans": 0.15,
                                   "freedom": -0.2,
                                   "warriors": 0.5}
        self.values["movement"] = {"travelers": 0.8,
                                   "sailors": 0.2,
                                   "traders": 0.55,
                                   "astrologists": 0.15,
                                   "builders": -0.4,
                                   "simplicity": 0.4,
                                   "materialists": -0.15,
                                   "freedom": 0.85,
                                   "naturalists": 0.2,
                                   "collectivists": -0.2,
                                   "warriors": 0.2}
        self.values["plantlife"] = {"agriculture": 0.25,
                                    "greed": 0.15,
                                    "naturalists": 0.25,
                                    "shamans": 0.05,
                                    "craftsmen": 0.4,
                                    "traders": 0.15,
                                    "builders": 0.4,
                                    "simplicity": -0.25,
                                    "collectivists": 0.15}
        self.values["nature"] = {"naturalists": 0.35,
                                 "shamans": 0.2,
                                 "agriculture": 0.25,
                                 "freedom": 0.25,
                                 "travelers": 0.15,
                                 "simplicity": 0.5,
                                 "collectivists": 0.15,
                                 "superstitious": 0.3,
                                 "astrologists": 0.1,
                                 "metallurgists": 0.4,
                                 "warriors": 0.05}
        self.values["growth"] = {"shamans": 0.1,
                                 "agriculture": 0.55,
                                 "naturalists": 0.45,
                                 "metallurgists": -0.2,
                                 "freedom": 0.15,
                                 "astrologists": -0.3,
                                 "collectivists": 0.35,
                                 "materialists": -0.1,
                                 "warriors": 0.05,
                                 "builders": 0.15}
        self.values["sky"] = {"travelers": 0.55,
                              "craftsmen": 0.3,
                              "traders": 0.25,
                              "superstitious": 0.5,
                              "metallurgists": 0.1,
                              "worshippers": 0.5,
                              "freedom": 0.55,
                              "shamans": 0.1,
                              "astrologists": 0.7,
                              "simplicity": 0.4,
                              "builders": 0.1}
        self.values["earth"] = {"metallurgists": 1.2,
                                "craftsmen": 0.9,
                                "traders": 0.45,
                                "materialists": 0.75,
                                "agriculture": 0.25,
                                "collectivists": 0.4,
                                "builders": 0.6,
                                "freedom": -0.25,
                                "worshippers": -0.15,
                                "greed": 0.6,
                                "superstitious": -0.25,
                                "warriors": 0.05,
                                "astrologists": -0.2}
        self.values["fields"] = {"agriculture": 0.75,
                                 "builders": 0.3,
                                 "materialists": 0.5,
                                 "naturalists": 0.35,
                                 "superstitious": -0.3,
                                 "simplicity": -0.3,
                                 "collectivists": 0.25,
                                 "warriors": 0.05}
        self.values["sunlight"] = {"worshippers": 0.9,
                                   "astrologists": 0.6,
                                   "naturalists": 0.15,
                                   "travelers": 0.25,
                                   "simplicity": 0.75,
                                   "freedom": 0.6,
                                   "materialists": -0.2,
                                   "warriors": 0.3,
                                   "builders": -0.1}
        self.values["ice"] = {"superstitious": 0.3,
                              "simplicity": -0.4,
                              "freedom": -0.4,
                              "travelers": -0.45,
                              "materialists": 0.4,
                              "sailors": 0.1,
                              "shamans": 0.45,
                              "greed": 0.7,
                              "metallurgists": 0.2,
                              "warriors": 0.4,
                              "agriculture": -0.2}
        self.values["fear"] = {"superstitious": 0.85,
                               "worshippers": 0.3,
                               "shamans": 0.8,
                               "freedom": -0.55,
                               "collectivists": 0.4,
                               "simplicity": -0.3,
                               "builders": 0.25,
                               "greed": 0.8,
                               "materialists": -0.25,
                               "warriors": 0.65}
        self.values["water"] = {"sailors": 2,
                                "simplicity": 0.3,
                                "freedom": 0.75,
                                "travelers": 0.6,
                                "builders": 0.15,
                                "craftsmen": 0.1,
                                "astrologists": 0.45,
                                "traders": 0.75,
                                "warriors": 0.1}

    def influences(self):
        self.influenceOutputs = {"sky": 0,
                                 "sunlight": 0,
                                 "fields": 0,
                                 "earth": 0,
                                 "ice": 0,
                                 "fear": 0,
                                 "water": 0,
                                 "growth": 0,
                                 "nature": 0,
                                 "plantlife": 0,
                                 "movement": 0,
                                 "darkness": 0,
                                 "swimming": 0,
                                 "food": 0}
        self.influences = {}
        self.influences["elevation"] = {"sky": 0.3,
                                        "sunlight": 0.15,
                                        "fields": 0.1,
                                        "earth": 0.1,
                                        "ice": 0.1,
                                        "fear": 0.05}
        self.influences["slope"] = {"sky": 0.1,
                                    "sunlight": 0.1,
                                    "earth": 0.1,
                                    "fear": 0.05}
        self.influences["temperature"] = {"sunlight": 0.2,
                                          "darkness": -0.1,
                                          "sky": 0.1,
                                          "movement": 0.2,
                                          "fear": 0.15,
                                          "water": 0.05,
                                          "food": 0.1}
        self.influences["vegetation"] = {"nature": 0.15,
                                         "food": 0.15,
                                         "darkness": 0.05,
                                         "plantlife": 0.2,
                                         "sunlight": -0.1,
                                         "fields": -0.1}
        self.influences["rainfall"] = {"water": 0.5,
                                       "growth": 0.2,
                                       "nature": 0.1,
                                       "sky": 0.25,
                                       "movement": 0.1,
                                       "darkness": 0.1,
                                       "plantlife": 0.1,
                                       "swimming": 0.15}
        self.influences["vegetation"] = {"plantlife": 0.6,
                                         "growth": 0.15,
                                         "nature": 0.35,
                                         "earth": 0.2,
                                         "fields": 0.15}
        self.influences["fertility"] = {"plantlife": 0.35,
                                        "growth": 0.25,
                                        "nature": 0.15,
                                        "earth": 0.3,
                                        "fields": 0.25}
        self.influences["metallicity"] = {"earth": 0.7,
                                          "darkness": 0.15,
                                          "nature": 0.1}
        self.influences["carnivores"] = {"darkness": 0.35,
                                         "food": 0.25,
                                         "fear": 0.4,
                                         "nature": 0.15}
        self.influences["herbivores"] = {"nature": 0.35,
                                         "food": 0.65,
                                         "growth": 0.3,
                                         "earth": 0.1,
                                         "plantlife": -0.1}
        self.influences["river"] = {"water": 2.5,
                                    "swimming": 1.0,
                                    "nature": 0.1,
                                    "food": 0.35,
                                    "growth": 0.1,
                                    "movement": 0.3}
        self.influences["water"] = {"water": 4.0,
                                    "swimming": 1.5,
                                    "sky": 0.15,
                                    "fear": 0.2,
                                    "darkness": 0.35,
                                    "food": 0.25,
                                    "movement": 0.15,
                                    "ice": 0.1}
        self.influences["forest"] = {"plantlife": 0.6,
                                     "darkness": 0.25,
                                     "nature": 0.4,
                                     "growth": 0.2,
                                     "earth": 0.1,
                                     "food": 0.3,
                                     "fields": 0.1}
        self.influences["desert"] = {"earth": 0.4,
                                     "sky": 0.3,
                                     "sunlight": 0.7,
                                     "food": -0.3,
                                     "fields": 0.35}
        self.influences["shrubland"] = {"food": 0.25,
                                        "fields": 0.4,
                                        "sky": 0.35,
                                        "sunlight": 0.2,
                                        "earth": 0.25,
                                        "nature": 0.1}
        self.influences["savanna"] = {"food": 0.2,
                                      "fields": 0.6,
                                      "sky": 0.4,
                                      "sunlight": 0.3,
                                      "nature": 0.2,
                                      "earth": 0.1}
        self.influences["tundra"] = {"food": -0.35,
                                     "fields": 0.35,
                                     "fear": 0.2,
                                     "darkness": 0.4,
                                     "sky": 0.2,
                                     "growth": -0.1,
                                     "earth": 0.3,
                                     "ice": 0.5}
        self.influences["mountain"] = {"food": -0.3,
                                       "fields": 0.1,
                                       "sky": 0.7,
                                       "sunlight": 0.35,
                                       "earth": 0.45,
                                       "fear": 0.2,
                                       "nature": 0.15,
                                       "ice": 0.5}
        self.influences["tropical"] = {"food": 0.35,
                                       "fear": 0.4,
                                       "darkness": 0.35,
                                       "fields": -0.25,
                                       "nature": 0.5,
                                       "plantlife": 0.6,
                                       "growth": 0.45,
                                       "earth": 0.1}
        self.influences["boreal"] = {"food": 0.1,
                                     "nature": 0.25,
                                     "growth": 0.1,
                                     "fields": 0.2,
                                     "darkness": 0.2,
                                     "plantlife": 0.25,
                                     "earth": 0.15,
                                     "ice": 0.15}
        self.influences["frozen"] = {"ice": 0.75,
                                     "earth": 0.2,
                                     "food": -0.3,
                                     "plantlife": -0.2,
                                     "fear": 0.2,
                                     "darkness": 0.2,
                                     "fields": 0.15,
                                     "water": 0.1}

    def placeCity(self, xx, yy, pop=50, culture=None, node=None):
        if node is not None:
            cityNode = self.nearestNode(xx, yy)
        else:
            cityNode = node
        if cityNode.biome != "water" and cityNode.city is None:
            newCity = City(cityNode, pop, culture, self)
            return 1
        else:
            return -1

    def randomCity(self):
        cityNode = self.atlas[math.floor(random.random() * len(self.atlas))]
        while (cityNode.biome == "water" or cityNode.city is not None or
               cityNode.x < 32 or cityNode.x > self.xDim - 32 or cityNode.y < 32 or cityNode.y > self.yDim - 32):
            if len(self.cities) >= 1:
                cityNode = self.atlas[math.floor(random.random() * len(self.atlas))]
                while cityNode.dist(self.nearestCity(cityNode.x, cityNode.y).node) < 32:
                    cityNode = self.atlas[math.floor(random.random() * len(self.atlas))]
            else:
                cityNode = self.atlas[math.floor(random.random() * len(self.atlas))]
        newCity = City(cityNode, pop=random.randint(12, 136), m=self)

    def scatterCities(self, n):
        for i in range(n):
            self.randomCity()

    def drawGraph(self, gui=None):
        visualAtlas = Image.new("HSV", (mapDimX, mapDimY), "white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawGraph(graphDraw)
        for n in self.atlas:
            n.drawPoint(graphDraw, 1, "red")
        visualAtlas.show()

    def drawElevation(self, pts=0, sea=1, gui=None):
        visualAtlas = Image.new("HSV", (mapDimX, mapDimY), "white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawElevation(graphDraw, self.sealevel * sea)
        for n in self.atlas:
            n.drawElevation(graphDraw)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw)
        visualAtlas.show()

    def drawLandmass(self, pts=0, nbrs=0, gui=None):
        visualAtlas = Image.new("HSV", (mapDimX, mapDimY), "white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawLandmass(graphDraw, self.sealevel)
        for n in self.atlas:
            n.drawLandmass(graphDraw, pts, nbrs)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawPath(graphDraw)
        visualAtlas.show()

    def drawWildlife(self, gui=None):
        visualAtlas = Image.new("HSV", (mapDimX, mapDimY), "white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawWildlife(graphDraw, self.sealevel)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw, self.xDim)
        visualAtlas.show()

    def displayNode(self, event):
        clickedNode = self.nearestNode(event.x, event.y)
        self.displayString.set(self.nodeInfo(clickedNode))
        self.displayNo = clickedNode
        cityNode = self.nearestCity(event.x, event.y)
        if cityNode.node.dist(Node(event.x, event.y)) < 8:
            self.displayString.set(self.nodeInfo(cityNode.node))
            self.displayNo = cityNode.node

    def redraw(self):
        visualAtlas = Image.new("HSV", (mapDimX, mapDimY), "white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawReal(graphDraw, self.sealevel)
        for n in self.atlas:
            n.drawReal(graphDraw, self.sealevel)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw, self.xDim)
        for c in self.cities:
            c.drawSelf(graphDraw)
        visualAtlas = visualAtlas.convert("RGB")
        visualAtlas.save(self.mapname, "GIF")
        photo = Image.open(self.mapname)
        self.img = ImageTk.PhotoImage(photo)
        self.lbl.configure(image=self.img)
        self.lbl.image = self.img
        if self.displayNo is not None:
            self.displayString.set(self.nodeInfo(self.displayNo))

    def updateResources(self):
        for r in self.resourceRegions:
            r.updateReg()

    def updateTerritory(self):
        for c in self.cities:
            c.node.culture = c.culture
            c.node.allegiance = 1 / c.population
        for p in self.atlas:
            p.updateAllegiance()
        for c in self.cities:
            c.node.culture = c.culture
            c.node.allegiance = 1 / c.population

    def updateDemogs(self):
        for c in self.cities:
            c.updateDemog()

    def nextTurn(self):
        self.updateResources()
        self.updateDemogs()
        self.updateTerritory()
        self.redraw()

    def cultureInfo(self):
        if self.displayNo is None:
            return -1
        if self.displayNo.culture is None:
            return -1
        self.displayCulture = self.displayNo.culture
        if self.infoGui is not None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        photo = Image.open(self.displayCulture.flag.filename)
        self.flagImg = ImageTk.PhotoImage(photo)
        self.flagLbl = Label(self.infoGui, image=self.flagImg)
        self.flagLbl.config(borderwidth=32)
        self.flagLbl.photo = photo
        self.flagLbl.pack()
        self.cultureString = StringVar()
        self.cultureString.set(self.displayCulture.cultureNotes())
        cdsc = Label(self.infoGui, textvariable=self.cultureString)
        cdsc.pack()

    def drawReal(self, gui=None):
        visualAtlas = Image.new("HSV", (mapDimX, mapDimY), "white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawReal(graphDraw, self.sealevel)
        for n in self.atlas:
            n.drawReal(graphDraw, self.sealevel)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw, self.xDim)
        for c in self.cities:
            c.drawSelf(graphDraw)
        if gui is None:
            visualAtlas = visualAtlas.convert("RGB")
            visualAtlas.save("map00.png", "PNG")
            visualAtlas.show()
        else:
            self.gui = gui
            self.displayString = StringVar()
            self.displayString.set("No node selected")
            self.infoScales()
            desc = Label(gui, textvariable=self.displayString)
            desc.pack(side=RIGHT, fill=Y)
            visualAtlas = visualAtlas.convert("RGB")
            self.mapname = self.cultures[0].language.genName() + ".gif"
            visualAtlas.save(self.mapname, "GIF")
            photo = Image.open(self.mapname)
            self.img = ImageTk.PhotoImage(photo)
            self.lbl = Label(gui, image=self.img)
            self.lbl.pack()
            self.lbl.bind("<Button-1>", self.displayNode)
            b0 = Button(gui, text="Next Turn", command=self.nextTurn)
            b0.pack(anchor=S, side=RIGHT)
            c1 = "orange red"
            b0.config(bg=c1, activebackground=c1, activeforeground=c1)
            b1 = Button(gui, text="Society Info", command=self.cultureInfo)
            b1.pack(anchor=S, side=RIGHT)
            c1 = "medium aquamarine"
            b1.config(bg=c1, activebackground=c1, activeforeground=c1)
            self.nextTurn()
            gui.mainloop()
