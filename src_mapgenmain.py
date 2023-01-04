# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 22:47:38 2017

@author: Bri

"""

import numpy as np
import random
import math
from tkinter import *
from PIL import Image, ImageDraw, ImageTk
from src_towngen import *
from src_events import *
from src_items import *
from src_tools import *
from src_facegen import *
from src_magic import *
import string

def lerp(t,a,b):
    return a + t * (b - a)

def smoothcurve(t):
    return t * t * (3. - 2. * t)

def lengthDirX(length, angle):
  radian_angle = math.radians(angle)
  return length * math.cos(radian_angle)

def lengthDirY(length, angle):
  radian_angle = math.radians(angle)
  return length * math.sin(radian_angle)

def A(dx, dy):
  return math.degrees( math.atan2(dy, dx) )

def drawCircle(drawer,x,y,radius,color):
    x1 = x-radius
    x2 = x+radius
    y1 = y-radius
    y2 = y+radius
    drawer.ellipse([(x1,y1),(x2,y2)],color)

def drawTrapezoid(drawer,x1,y1,x2,y2,r1,r2,color):
    directAngle = A(x2-x1,y2-y1)
    pAngle = directAngle-90
    pAngle2 = pAngle-180
    p1 = (x1+lengthDirX(r1,pAngle),y1+lengthDirY(r1,pAngle))
    p2 = (x1+lengthDirX(r1,pAngle2),y1+lengthDirY(r1,pAngle2))
    p3 = (x2+lengthDirX(r2,pAngle2),y2+lengthDirY(r2,pAngle2))
    p4 = (x2+lengthDirX(r2,pAngle),y2+lengthDirY(r2,pAngle))
    drawer.polygon([p1,p2,p3,p4],color,color)
    
def centroidnp(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    return sum_x/length, sum_y/length

def unitCos(x):
    return (math.cos(x*math.pi)+1)/2

def distMod(x,maxDist):
    return unitCos(x/maxDist)

def stick(val,minimum,maximum):
    if abs(val-minimum) < abs(val-maximum):
        return minimum
    else:
        return maximum

def strDivider(length):
    n = ""
    for l in range(length):
        n += "_"
    return n
            
from scipy.spatial import Voronoi
def relaxLloyd(pts,strength):
    for i in range(strength):
        vor = Voronoi(pts)
        newpts = []
        for idx in range(len(vor.points)):
            pt = vor.points[idx,:]
            region = vor.regions[vor.point_region[idx]]
            if -1 in region:
                newpts.append(pt)
            else:
                vxs = np.asarray([vor.vertices[i,:] for i in region])
                newpt = centroidnp(vxs)
                newpts.append(newpt)
        pts = np.array(newpts) 
    
class Node:
    def __init__(self,xx,yy,m=None):
        self.x = xx
        self.y = yy
        self.neighbors = []
        self.landmass = None
        self.bodyWater = None
        self.waterdistance = 0
        self.river = None
        self.region = None
        self.city = None
        self.culture = None
        self.allegiance = 0
        self.resourceRegion = None
        self.resourceDist = 0
        self.key = 0
        self.roads = []
        self.entities = []
        self.myMap=m
        self.name = str(math.floor(xx+(yy*437)))
    def coords(self):
        tupleVert = (self.x,self.y)
        return tupleVert
    def dist(self,n):
        distX = abs(self.x-n.x)
        distY = abs(self.y-n.y)
        dist = math.sqrt((distX**2)+(distY**2))
        return dist
    def isLinked(self,nNode):
        if nNode in self.neighbors:
            return 1
        else:
            return 0
    def link(self,newNeighbor):
        if self.dist(newNeighbor) < 100 or len(self.neighbors) == 0:
            self.neighbors.append(newNeighbor)
            newNeighbor.neighbors.append(self)
    def linkRoads(self,newNeighbor):
        if newNeighbor not in self.neighbors:
            return 0
        if newNeighbor not in self.roads:
            self.roads.append(newNeighbor)
        if self not in newNeighbor.roads:
            newNeighbor.roads.append(self)
    def getKey(self):
        return self.key
    def watery(self,sealevel=0):
        if self.elevation < sealevel or self.bodyWater != None:
            return 1
        else:
            return 0
    def hasWaterNeighbor(self,sealevel):
        for n in self.neighbors:
            if n.watery(sealevel) == 1:
                if n.x > 0 and n.x < 960 and n.y > 0 and n.y < 960:
                    return 1
        return 0
    def nearestNeighbor(self,n):
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
        rise = maximum.elevation-minimum.elevation
        run = self.dist(minimum)+self.dist(maximum)
        if run != 0:
            m = rise/run
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
        self.rSlope = self.slope*self.slopeDir
        return self.rSlope
    def riverNext(self):
        if self.river == None:
            return None
        selfIndex = self.river.nodes.index(self)
        if selfIndex == len(self.river.nodes)-1:
            return None
        else:
            return self.river.nodes[selfIndex+1]
    def riverPrevious(self):
        if self.river == None:
            return None
        selfIndex = self.river.nodes.index(self)
        if selfIndex == 0:
            return None
        else:
            return self.river.nodes[selfIndex-1]
    def waterdist(self,sealevel):
        self.getSlope()
        if self.watery() == 1 or self.river != None:
            self.waterdistance = 0
        elif self.hasWaterNeighbor(sealevel) == 1:
            self.waterdistance = 1
        else:
            dd = [n.waterdistance for n in self.neighbors]
            self.waterdistance = min(dd) + 1 + (random.random()*0.4) + clamp(self.slope*100,0,1)
        if self.x < 0 or self.y < 0 or self.x > 960 or self.y > 960:
            self.waterdistance = 1000
    def smooth(self):
        nbrs = []
        nbrs.append(self.elevation)
        for i in self.neighbors:
            nbrs.append(i.elevation)
        self.elevation = sum(nbrs)/len(nbrs)
    def defaultRoads(self):
        self.roads = []
    def setVegetation(self):
        tempFitness = 1-abs(0.5-self.temp)
        elevationFitness = 1-abs(0.45-self.elevation)
        fertilityFitness = self.fertility-(self.metallicity*0.4)
        rainFitness = 1-abs(self.rainfall-0.8)
        vegFitness = ((tempFitness)*(elevationFitness+0.5)*(fertilityFitness+0.5)*(rainFitness))*45
        self.vegetation = clamp(vegFitness/128,0,1)
    def setBiome(self,sl):
        if self.elevation > 0.864:
                self.biome = "mountain"
        elif self.temp < 0.18:
            if self.rainfall < 0.2:
                self.biome = "frost"
            elif self.rainfall < 0.36:
                self.biome = "frost"
            elif self.rainfall < 0.55:
                self.biome = "shrubland"
            else:
                self.biome = "boreal forest"
        elif self.temp < 0.34:
            if self.rainfall < 0.1:
                self.biome = "tundra"
            elif self.rainfall < 0.16:
                self.biome = "shrubland"
            elif self.rainfall < 0.5:
                self.biome = "boreal forest"
            else:
                self.biome = "forest"
        elif self.temp < 0.58:
            if self.rainfall < 0.03:
                self.biome = "desert"
            elif self.rainfall < 0.06:
                self.biome = "savanna"
            elif self.rainfall < 0.1:
                self.biome = "shrubland"
            elif self.rainfall < 0.6:
                self.biome = "forest"
            else:
                self.biome = "tropical forest"
        else:
            if self.rainfall < 0.04:
                self.biome = "desert"
            elif self.rainfall < 0.1:
                self.biome = "savanna"
            else:
                self.biome = "tropical forest"
        if self.watery() == 1:
            self.biome = "water"
    def claim(self,n,sealevel=0.4):
        if n.resourceRegion != None and n.resourceRegion.culture != self.culture:
            return -1
        n.culture = self.culture
        inc = 1
        if n.watery() == 1:
            inc *= 16
        if self.dist(n) > 32:
            inc *= 16
        n.allegiance = self.allegiance+inc
        if self.culture.name not in n.region.culturalNames:
            n.region.culturalNames[self.culture.name] = self.culture.language.genName()
        if n.landmass != None:
            if self.culture.name not in n.landmass.culturalNames:
                n.landmass.culturalNames[self.culture.name] = self.culture.language.genName()
        if n.river != None:
            if self.culture.name not in n.river.culturalNames:
                n.river.culturalNames[self.culture.name] = self.culture.language.genName()
    def updateAllegiance(self,sealevel=0.4):
        if self.culture != None:
            chance = 1/clamp(self.allegiance,0.00001,512)
            for n in self.neighbors:
                if n.culture != None and len(self.culture.cultureOpinions) > 0 and n.culture != self.culture:
                    self.culture.cultureOpinions[n.culture.name].knowledge = 2
                roll = random.random()
                if roll <= chance:
                    self.claim(n,sealevel)
    def structure(self):
        structureSeed = int(self.name)
        possibleStructures = ["farm","inn","brothel","factory","mine","fort","farm","farm","farm","farm","farm","inn","mine"]
        if self.city != None:
            return None
        if self.resourceRegion == None:
            return None
        if self.watery() == 1:
            return None
        if structureSeed % 2 == 0:
            return None
        return possibleStructures[structureSeed % len(possibleStructures)]
    def toString(self):
        self.selfString = "("
        self.selfString += str(self.x)
        self.selfString += ", "
        self.selfString += str(self.y)
        self.selfString += ")"
        return self.selfString
    def printSelf(self):
        print(self.toString())
    def drawPoint(self,drawer,radius,color):
        drawCircle(drawer,self.x,self.y,radius,color)
    def drawElevation(self,drawer,pts=0):
        col = math.floor(self.elevation*255)
        dCol = (0,0,col)
        if self.elevation >= 1:
            dCol = "red"
        if self.landmass != None:
            for n in self.neighbors:
                if n.landmass == self.landmass:
                    drawer.line([self.coords(),n.coords()],dCol,3)
        if pts == 1:
            self.drawPoint(drawer,1,dCol)
    def drawLandmass(self,drawer,pts=0,drawNeighbors=0):
        if self.landmass == None:
            dCol = "black"
        else:
            dCol = self.landmass.color
        if self.landmass != None:
            for n in self.neighbors:
                if n.landmass == self.landmass:
                    drawer.line([self.coords(),n.coords()],dCol,3)
        if self.landmass != None and drawNeighbors == 1:
            for n in self.neighbors:
                drawer.line([self.coords(),n.coords()],dCol,1)
        if pts == 1:
            self.drawPoint(drawer,1,"black")
    def drawReal(self,drawer,sl):
        dCol = self.biomeColor
        if self.landmass != None and 0:
            for n in self.neighbors:
                if n.landmass == self.landmass:
                    drawer.line([self.coords(),n.coords()],dCol,3)
    def drawRoads(self,drawer,roadCol):
        for n in self.neighbors:
            if n in self.roads:
                drawer.line([self.coords(),n.coords()],roadCol,2)
    def drawTownGen(self):
        self.townGen = Town(self,self.myMap,self.name)
        townImg = Image.new("HSV",(self.townGen.xDim,self.townGen.yDim),(0,0,255))
        graphDraw = ImageDraw.Draw(townImg)
        self.townGen.drawSelf(graphDraw)
        self.townGen.drawRoads(townImg)
        townImg = townImg.convert("RGB")
        fileName = ("./generated/"+self.townGen.mapName)
        townImg.save(self.townGen.mapName,"GIF")
        

class Triangle:
    def __init__(self,aa,bb,cc,m):
        self.verts = [None,None,None]
        self.verts[0] = aa
        self.verts[1] = bb
        self.verts[2] = cc
        self.neighbors = []
        self.myMap = m
    def sharesNeighbors(self,other):
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
    def drawGraph(self,drawer):
        drawer.line([self.verts[0].coords(),self.verts[1].coords()],"black",1)
        drawer.line([self.verts[1].coords(),self.verts[2].coords()],"black",1)
        drawer.line([self.verts[2].coords(),self.verts[0].coords()],"black",1)
    def drawElevation(self,drawer,sealevel=0):
        elevationList = [self.verts[f].elevation for f in range(len(self.verts))]
        elevation = sum(elevationList)/3
        col = math.floor(elevation*255)
        dCol = (0,0,col)
        underwater = 0
        for f in self.verts:
            if f.watery() == 1:
                underwater = 1
        if underwater == 1:
            dCol = (142,128,clamp(col,64,76))
        drawer.polygon([self.verts[0].coords(),self.verts[1].coords(),self.verts[2].coords()],fill=dCol,outline=dCol)
    def drawWildlife(self,drawer,sealevel=0):
        carnivores = math.floor(128*sum([self.verts[f].carnivores for f in range(len(self.verts))])/3)
        herbivores = math.floor(128*sum([self.verts[f].herbivores for f in range(len(self.verts))])/3)
        dCol = (128+(carnivores)-(herbivores),(carnivores+herbivores)*2,128)
        underwater = 0
        for f in self.verts:
            if f.watery() == 1:
                underwater = 1
        if underwater == 1:
            dCol = (142,128,70)
        drawer.polygon([self.verts[0].coords(),self.verts[1].coords(),self.verts[2].coords()],fill=dCol,outline=dCol)
    def drawLandmass(self,drawer,sealevel=0):
        landmass = None
        for p in self.verts:
            if p.landmass != None:
                landmass = p.landmass
        if landmass != None:
            dCol = landmass.color
        else:
            dCol = "black"
        underwater = 0
        for f in self.verts:
            if f.watery() == 1:
                underwater = 1
        if underwater == 1:
            dCol = (142,64,64)
        drawer.polygon([self.verts[0].coords(),self.verts[1].coords(),self.verts[2].coords()],fill=dCol,outline=dCol)
    def drawTerritory(self,drawer,sealevel=0):
        t = [x.culture for x in self.verts]
        c = max(set(t),key=t.count)
        h = 0
        s = 0
        v = 128
        if c != None:
            h = c.bannerColor[0]
            s = clamp(c.bannerColor[1],32,255)
            v = clamp(c.bannerColor[2],96,192)
        underwater = 0
        for f in self.verts:
            if f.watery() == 1:
                underwater = 1
        if underwater == 1:
            v = 32
        dCol = (h,s,v)
        tt = [x.resourceRegion for x in self.verts]
        cc = max(set(tt),key=tt.count)
        if cc != None:
            dCol = (0,0,0)
        drawer.polygon([self.verts[0].coords(),self.verts[1].coords(),self.verts[2].coords()],fill=dCol,outline=dCol)
    def drawReal(self,drawer,sl):
        elevationList = [self.verts[f].elevation for f in range(len(self.verts))]
        elevation = sum(elevationList)/3
        col = math.floor(((elevation*255)+64)/2)
        underwater = 0
        avgHue = math.floor(sum([self.verts[f].biomeColor[0] for f in range(len(self.verts))])/3)
        avgSat = math.floor(sum([self.verts[f].biomeColor[1] for f in range(len(self.verts))])/3)
        avgVal = math.floor(sum([self.verts[f].biomeColor[2] for f in range(len(self.verts))])/3)
        dCol = (avgHue,avgSat,avgVal)
        for f in self.verts:
            if f.watery() == 1:
                underwater = 1
        if underwater == 1:
            v = clamp(col,self.myMap.waterValueMin,self.myMap.waterValueMax)
            dCol = (self.myMap.waterHue,self.myMap.waterSaturation,v)
        drawer.polygon([self.verts[0].coords(),self.verts[1].coords(),self.verts[2].coords()],fill=dCol,outline=dCol)

class River:
    def __init__(self,root,length,landms):
        self.nodes = []
        self.addNode(root)
        self.landmass = landms
        current = root
        j = 0
        while j < length:
            choice = current
            m = 1000
            for k in current.neighbors:
                kr = 0
                if k.river != None and k.river != self:
                    kr = 1
                for q in k.neighbors:
                    if q.river != None and q.river != self:
                        kr = 1
                if (k.elevation < m and k.elevation >= current.elevation and kr == 0 and not k.hasWaterNeighbor(self.landmass.sealevel)):
                    choice = k
                    m = k.elevation
            if choice == current:
                j = length
            else:
                self.addNode(choice)
            j+=1
            current = choice
        self.culturalNames = {}
    def removeRiver(self):
        for n in self.nodes:
            n.river = None
            self.nodes.remove(n)
    def addNode(self,newNode):
        self.nodes.append(newNode)
        newNode.river = self
    def drawPath(self,drawer):
        dCol = (142,64,64)
        drawNodes = []
        for n in self.nodes:
            drawNodes.append(n.coords())
        drawer.line(drawNodes,dCol,2)
    def drawRiver(self,drawer,xDim):
        for l in self.nodes:
            l.getSlope()
        dCol = (self.landmass.myMap.waterHue,self.landmass.myMap.waterSaturation,self.landmass.myMap.waterValueMax)
        for i in range(len(self.nodes)-1):
            dCol = (self.landmass.myMap.waterHue,self.landmass.myMap.waterSaturation,self.landmass.myMap.waterValueMax)
            n = self.nodes[i]
            n1 = self.nodes[i+1]
            scale = xDim/2
            w = clamp((1/n.slope)/scale,0.5,2)
            w1 = clamp((1/n1.slope)/scale,0.5,2)
            drawCircle(drawer,n.x,n.y,w,dCol)
            drawCircle(drawer,n1.x,n1.y,w1,dCol)
            drawTrapezoid(drawer,n.x,n.y,n1.x,n1.y,w,w1,dCol)

class bodyWater:
    def __init__(self,rootNode,sLevel,maxsize=100000):
        self.sealevel = sLevel
        self.root = rootNode
        self.nodes = []
        self.addNode(rootNode)
        self.fill()
        self.culturalNames = {}
        self.maxsize = maxsize
    def addNode(self,p):
        if p not in self.nodes:
            self.nodes.append(p)
            p.bodyWater = self
    def fill(self):
        for p in self.nodes:
            for k in p.neighbors:
                while k not in self.nodes and k.elevation < self.sealevel:
                    self.addNode(k)

class Landmass:
    def __init__(self,rootNode,m):
        self.myMap = m
        self.sealevel = self.myMap.sealevel
        self.color = self.landmassColor()
        self.root = rootNode
        self.nodes = []
        if rootNode.watery() == 1:
            return 0
        self.addNode(rootNode)
        self.boundary = []
        self.rivers = []
        self.fill()
        self.size = len(self.nodes)
        self.centermass()
        self.culturalNames = {}
        self.landmassType = self.lType(self.size)
    def centermass(self):
        xTotal = sum([p.x for p in self.nodes])
        yTotal = sum([p.y for p in self.nodes])
        xx = xTotal/self.size
        yy = yTotal/self.size
        self.centroid = Node(xx,yy,self.myMap)
    def lType(self,s):
        if s <= 6:
            t = "islet"
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
        h = math.floor(random.random()*255)
        s = 128+math.floor(random.random()*128)
        v = 128+math.floor(random.random()*128)
        return (h,s,v)
    def addNode(self,p):
        if p not in self.nodes:
            self.nodes.append(p)
            p.landmass = self
    def addBoundary(self,p):
        if p not in self.boundary:
            self.boundary.append(p)
    def fill(self):
        for p in self.nodes:
            for k in p.neighbors:
                if k.watery() == 1:
                    self.addBoundary(p)
                while k not in self.nodes and k.watery() == 0:
                    self.addNode(k)
    def addRiver(self,length):
        root = random.choice(self.boundary)
        riverLen = random.randint(length/4,length)
        newRiver = River(root,riverLen,self)
        self.rivers.append(newRiver)
    def cullStreams(self,minLen):
        for r in self.rivers:
            if len(r.nodes) < minLen:
                self.rivers.remove(r)
                r.removeRiver()

class Region:
    def __init__(self,rootNode):
        self.root = rootNode
        self.nodes = []
        self.culturalNames = {}
        self.addNode(self.root)
        self.biome = self.root.biome
        self.fill()
        if self.biome == "boreal" or self.biome == "tropical":
            self.biome += " forest"
        if self.biome == "frost":
            self.biome += " tundra"
        if self.biome == "water":
            bodysize = len(self.nodes)
            if bodysize > 2048:
                self.biome = "ocean"
            elif bodysize > 512:
                self.biome = "sea"
            else:
                self.biome = "lake"
    def addNode(self,p):
        if p not in self.nodes:
            self.nodes.append(p)
            p.region = self
    def fill(self):
        for p in self.nodes:
            for k in p.neighbors:
                if k.biome == self.biome:
                    self.addNode(k)

class Influence:
    def __init__(self,mMap,myNode,root):
        self.node = myNode
        self.rootNode = root
        self.myMap = mMap
        self.influenceOutput = {}
        self.envInfluences = {}
        for o in self.myMap.influenceOutputs.keys():
            self.influenceOutput[o] = 0
        self.translateInfluence()
    def maxInf(self):
        for k in self.influenceOutput:
            maxkey = k
        for k in self.influenceOutput:
            if self.influenceOutput[k] > self.influenceOutput[maxkey]:
                maxkey = k
        return maxkey
    def influencesMain(self,n):
        mInf = {}
        for f in range(n):
            maxkey = self.maxInf()
            mInf[maxkey] = self.influenceOutput[maxkey]
            del self.influenceOutput[maxkey]
        for i in mInf.keys():
            self.influenceOutput[i] = self.mInf[i]
        return mInf
    def setOutput(self,influenceType,strength):
        if influenceType in self.envInfluences.keys():
            self.envInfluences[influenceType] += strength
        else:
            self.envInfluences[influenceType] = strength
        for o in self.myMap.influenceOutputs:
            if o in self.myMap.influences[influenceType]:
                output = strength*self.myMap.influences[influenceType][o]
                self.influenceOutput[o] += output
    def translateInfluence(self):
        if self.rootNode == 1:
            for p in self.node.neighbors:
                ni = Influence(self.myMap,p,0)
                for j in ni.influenceOutput.keys():
                    self.influenceOutput[j] += ni.influenceOutput[j]/len(self.node.neighbors)
        strength = 1
        inf = self.node.biome
        self.setOutput(inf,strength)
        strength = self.node.herbivores*10
        inf = "carnivores"
        self.setOutput(inf,strength)
        strength = self.node.carnivores*10
        inf = "herbivores"
        self.setOutput(inf,strength)
        strength = self.node.temp
        inf = "temperature"
        self.setOutput(inf,strength)
        strength = self.node.elevation
        inf = "elevation"
        self.setOutput(inf,strength)
        strength = self.node.slope*100
        inf = "hills"
        self.setOutput(inf,strength)
        strength = self.node.vegetation
        inf = "vegetation"
        self.setOutput(inf,strength)
        strength = self.node.metallicity
        inf = "metallicity"
        self.setOutput(inf,strength)
        strength = self.node.fertility
        inf = "fertility"
        self.setOutput(inf,strength)
        strength = (self.node.dist(self.myMap.north)/self.myMap.xDim)
        inf = "latitude"
        self.setOutput(inf,strength)
        strength = self.node.rainfall
        inf = "rainfall"
        self.setOutput(inf,strength)
        
class Values:
    def __init__(self,m,inf):
        self.myMap = m
        self.influences = inf
        self.valuesOutput = {}
        for v in self.myMap.valuesOutputs.keys():
            self.valuesOutput[v] = 0
        self.translateValues()
        self.mainValues = self.valuesMain(5)
    def maxval(self):
        for k in self.valuesOutput.keys():
            maxkey = k
        for k in self.valuesOutput.keys():
            if self.valuesOutput[k] > self.valuesOutput[maxkey]:
                maxkey = k
        return maxkey
    def valuesMain(self,n):
        mVals = {}
        for f in range(n):
            maxkey = self.maxval()
            mVals[maxkey] = self.valuesOutput[maxkey]
            del self.valuesOutput[maxkey]
        for v in mVals.keys():
            self.valuesOutput[v] = mVals[v]
        return mVals
    def translateValues(self):
        for q in self.influences.influenceOutput.keys():
            modifier = self.influences.influenceOutput[q]
            roll = random.random()
            if roll > 0.97:
                modifier = 8
            elif roll < 0.04:
                modifier = 0.01
            for v in self.myMap.values[q].keys():
                self.valuesOutput[v] += self.myMap.values[q][v]*modifier*random.uniform(0.7,1.4)

class City:
    def __init__(self,n,pop=50,cltr=None,m=None):
        self.myMap = m
        self.myMap.cities.append(self)
        self.node = n
        self.node.city = self
        self.population = pop
        self.popThresholds = [40,200,1000,10000,50000]
        self.culture = None
        if cltr == None:
            self.culture = Culture(self.node,m)
        else:
            self.culture = cltr
        self.name = self.culture.language.genName()
        self.cityType = self.cType(self.population)
        self.region = self.node.region
        self.node.culture = self.culture
        self.node.allegiance = 1/self.population
        if self.culture.name not in self.node.region.culturalNames:
            self.node.region.culturalNames[self.culture.name] = self.culture.language.genName()
        if self.culture.name not in self.node.landmass.culturalNames:
            self.node.landmass.culturalNames[self.culture.name] = self.culture.language.genName()
        for q in self.node.neighbors:
            if q.region.biome == "ocean" or q.region.biome == "lake" or q.region.biome == "sea":
                if self.culture.name not in q.region.culturalNames:
                    q.region.culturalNames[self.culture.name] = self.culture.language.genName()
        if self.node.river != None:
            if self.culture.name not in self.node.river.culturalNames:
                self.node.river.culturalNames[self.culture.name] = self.culture.language.genName()
        reg = ResourceRegion(self,self.myMap)
        self.rawResources = [0,0]
        self.industrialProduction = 1
        self.foodProduction = 1
        self.age = 0
        self.roads = []
        self.kind = "city"
        e = Event(self.culture.myMap,a=self.age,kind="founding",sub=self,actrs=[self.culture.leader])
        e.importance = random.randint(10,40)+clamp(math.floor(math.sqrt(self.population)),0,25)
    def distanceToCity(self,other):
        return self.node.dist(other.node)
    def happiness(self):
        h = 1
        h = h-(math.log10(self.population)/100)
        if self.culture.society in ["Hegemony","Empire","Imperium"]:
            h = h-0.1
        h = h-abs(self.node.temp-0.4)**2
        m = self.culture.value.mainValues
        if "simplicity" in m:
            h += 0.1
        if "freedom" in m:
            h += 0.05
        return h
    def updateDemog(self):
        self.age += 1
        q = None
        for p in self.node.entities:
            if p.kind == "army" and p.profession == "guard infantry":
                q = p
        garrmod = 1
        if self.culture.society in ["Raiders","Pirates"]:
            garrmod = 1.6
        if self.culture.society in ["Shamanistic warriors","Religious zealots"]:
            garrmod = 1.4
        if self.culture.society in ["Hegemony","Empire","Imperium","Nation-state"]:
            garrmod = 1.2
        garrtarget = math.floor(clamp(self.population*0.08*garrmod,1,1+(2*(self.industrialProduction+self.foodProduction))))
        if q == None:
            q = Population(c=self.culture,n=None,t="",a=None,p=garrtarget,kind="army",node=self.node,prf="guard infantry")
            self.population = clamp(self.population-garrtarget,1,10000000)
        if q.number != garrtarget:
            d = garrtarget-q.number
            q.number = clamp(q.number + d,1,10000000)
            self.population = clamp(self.population-d,1,10000000)
        rscShare = self.population/self.node.resourceRegion.totalPop
        rscMax = [0,0]
        rscMax[0] = rscShare*self.node.resourceRegion.resources[0]
        rscMax[1] = rscShare*self.node.resourceRegion.resources[1]
        mpo = random.uniform(0.097,0.101)   # Maximum personal output (max resources production per person)
        m = self.culture.value.mainValues
        if "agriculture" in m:
            mpo *= 1.07
        if "simplicity" in m:
            mpo *= 0.85
        if "warriors" in m:
            mpo *= 1.04
        if "builders" in m:
            mpo *= 1.04
        if "metallurgists" in m:
            mpo *= 1.04
        if "craftsmen" in m:
            mpo *= 1.04
        if "collectivists" in m:
            mpo *= 1.03
        mpo = mpo*math.sqrt(self.culture.tech["production"])
        self.foodProduction = min(self.population*mpo,rscMax[0])
        self.industrialProduction = min(self.population*mpo,rscMax[1])
        mpc = random.uniform(0.079,0.083)     # Maximum personal consumption (max food needed per person)
        mpc -= 0.002*(math.log10(clamp(self.population,1,1000000))+2)
        # As the population grows, need less food per person due to economies of scale
        mpc += 0.002*(math.log2(clamp(len(self.node.resourceRegion.nodes),1,1000000))-2)
        # As the resource region grows, need more food per person due to transporation distance
        if "collectivists" in m:
            mpc -= 0.001
        if "freedom" in m:
            mpc += 0.001
        if "simplicity" in m:
            mpc -= 0.001
        if "greed" in m:
            mpc += 0.001
        mpc = mpc/math.sqrt(self.culture.tech["government"])
        self.foodConsumption = mpc*self.population
        diff = self.foodProduction-self.foodConsumption
        growth = clamp(diff/mpc,-self.population*0.12,self.population*0.06)
        if self.population < 20:
            growth = clamp(growth,1,100)
        self.population = math.ceil(self.population*0.99)  # Age related death
        self.population = clamp(math.floor(self.population+growth+random.choice([-1,0,0,0,0,1])),1,10000000)
        self.cityType = self.cType(self.population)
        self.diaspora()
    def diaspora(self):
        self.threshold = 1500*(0.99**self.age)
        roll = random.random()
        minRoll = 0.7
        superRoll = 0.98
        m = self.culture.value.mainValues
        if "simplicity" in m:
            self.threshold *= 1.15
            roll = roll/2
        if "travelers" in m:
            self.threshold *= 0.8
            roll = (roll+1)/2
        if "freedom" in m:
            self.threshold *= 0.9
        if "naturalists" in m:
            self.threshold *= 0.9
        if "materialists" in m:
            self.threshold *= 1.1
        if "agriculture" in m:
            self.threshold *= 1.1
        if "collectivists" in m:
            self.threshold *= 1.15
        if "builders" in m:
            self.threshold *= 1.2
        if "warriors" in m:
            self.threshold *= 1.2
        if (((self.population > self.threshold and roll > minRoll) or roll > superRoll)
        and self.age > 12 and (random.random() < 1/(len(self.culture.cultureCities())))):
            self.migrate()
    def migrate(self):
        emigrants = math.ceil(self.population*random.uniform(0.1,0.5))
        rng = 96
        xx = clamp(self.node.x + random.randint(-rng,rng),16,self.myMap.xDim)
        yy = clamp(self.node.y + random.randint(-rng,rng),16,self.myMap.yDim)
        n = self.myMap.nearestNode(xx,yy)
        if ((n.culture == self.culture or n.culture == None) 
            and n.city == None and n.resourceRegion == None and n.landmass != None):
            cc = City(n,pop=emigrants,cltr=self.culture,m=self.myMap)
            bcount = random.randint(19,23)
            builders = Population(c=self.culture,t="builders",a=random.randint(22,32),p=bcount,kind="group",node=self.node,prf="roadbuilder")
            builders.setPath(n)
            emigrants += bcount
            self.population = clamp(self.population-emigrants,1,10000000)
            self.age = 0
    def cType(self,p):
        if p <= self.popThresholds[0]:
            c = synonym("camp",seedNum(self.name))
        elif p <= self.popThresholds[1]:
            c = synonym("village",seedNum(self.name))
        elif p <= self.popThresholds[2]:
            c = synonym("township",seedNum(self.name))
        elif p <= self.popThresholds[3]:
            c = "town"
        elif p <= self.popThresholds[4]:
            c = "city"
        else:
            c = "metropolis"
        return c.capitalize()
    def governanceName(self):
        if self.culture.origin == self.node:
            return self.culture.nameOfCapital()
        else:
            return self.culture.nameOfTownHall(self.name)
    def justName(self):
        return self.name
    def nameFull(self):
        return self.cType(self.population) + " " + self.name
    def cityInfo(self):
        n = self.name + " (" + self.cityType + ")\n"
        n += "Population: " + str(self.population)
        return n
    def cultureInfo(self):
        info = self.culture.information()
        return info
    def drawTent(self,drawer,xx,yy,col,out):
        p0 = (xx,yy-6)
        p1 = (xx-4,yy+2)
        p2 = (xx+4,yy+2)
        drawer.polygon([p0,p1,p2],outline=out,fill=col)
        p0 = (xx,yy-1)
        p1 = (xx-1,yy+2)
        p2 = (xx+1,yy+2)
        drawer.polygon([p0,p1,p2],outline=out,fill=out)
        p0 = (xx,yy-6)
        p1 = (xx-2,yy-8)
        p2 = (xx+2,yy-8)
        drawer.line([p0,p1],fill=out,width=1)
        drawer.line([p0,p2],fill=out,width=1)
    def drawHut(self,drawer,xx,yy,col,out):
        p0 = (xx+3,yy+2)
        p1 = (xx-3,yy+2)
        p2 = (xx-3,yy-3)
        p3 = (xx,yy-6)
        p4 = (xx+3,yy-3)
        drawer.polygon([p0,p1,p2,p3,p4],outline=out,fill=col)
        p1 = (xx-1,yy-1)
        p2 = (xx+1,yy+2)
        drawer.rectangle([p1,p2],outline=out,fill=out)
        p0 = (xx,yy-6)
        p1 = (xx+5,yy-1)
        p2 = (xx-5,yy-1)
        drawer.line([p0,p1],fill=out,width=1)
        drawer.line([p0,p2],fill=out,width=1)
    def drawVillage(self,drawer,xx,yy,col,out):
        p0 = (xx,yy-5)
        p1 = (xx-2,yy-7)
        p2 = (xx-4,yy-5)
        p3 = (xx-4,yy+3)
        p4 = (xx,yy+3)
        drawer.polygon([p0,p1,p2,p3,p4],outline=out,fill=col)
        p0 = (xx-6,yy-3)
        p2 = (xx+2,yy-3)
        drawer.line([p1,p0],fill=out,width=1)
        drawer.line([p1,p2],fill=out,width=1)
        p3 = (xx-2,yy-4)
        p4 = (xx-2,yy-3)
        drawer.line([p4,p3],fill=out,width=1)
        p0 = (xx,yy-1)
        p1 = (xx+2,yy-3)
        p2 = (xx+4,yy-1)
        p3 = (xx+4,yy+3)
        p4 = (xx,yy+3)
        drawer.polygon([p0,p1,p2,p3,p4],outline=out,fill=col)
        p0 = (xx-2,yy+1)
        p2 = (xx+6,yy+1)
        drawer.line([p1,p0],fill=out,width=1)
        drawer.line([p1,p2],fill=out,width=1)
        p3 = (xx+2,yy+3)
        p4 = (xx+2,yy+1)
        drawer.line([p4,p3],fill=out,width=1)
    def drawTown(self,drawer,xx,yy,col,out):
        p0 = (xx-4,yy+3)
        p1 = (xx-4,yy-2)
        p2 = (xx,yy-5)
        p3 = (xx+4,yy-2)
        p4 = (xx+4,yy+3)
        drawer.polygon([p0,p1,p3,p4],outline=out,fill=col)
        drawer.polygon([p1,p2,p3],outline=out,fill=out)
        p0 = (xx-2,yy+4)
        p1 = (xx-2,yy-5)
        p2 = (xx,yy-7)
        p3 = (xx+2,yy-5)
        p4 = (xx+2,yy+4)
        drawer.polygon([p0,p1,p3,p4],outline=out,fill=col)
        drawer.polygon([p1,p2,p3],outline=out,fill=out)
        p0 = (xx,yy)
        p1 = (xx,yy+3)
        drawer.line([p1,p0],fill=out,width=1)
    def drawCity(self,drawer,xx,yy,col,out):
        p0 = (xx-2,yy-6)
        p1 = (xx-2,yy+2)
        p2 = (xx+2,yy+2)
        p3 = (xx+2,yy-6)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx-5,yy+1)
        p1 = (xx-5,yy-1)
        p2 = (xx+5,yy-1)
        p3 = (xx+5,yy+1)
        p4 = (xx,yy+4)
        drawer.polygon([p0,p1,p2,p3,p4],outline=out,fill=col)
        p0 = (xx+6,yy-3)
        p1 = (xx+6,yy+2)
        p2 = (xx+4,yy+2)
        p3 = (xx+4,yy-3)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx-6,yy-3)
        p1 = (xx-6,yy+2)
        p2 = (xx-4,yy+2)
        p3 = (xx-4,yy-3)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx,yy-1)
        p1 = (xx,yy+3)
        drawer.line([p1,p0],fill=out,width=1)
    def drawMetropolis(self,drawer,xx,yy,col,out):
        p0 = (xx-5,yy-1)
        p1 = (xx-6,yy+2)
        p2 = (xx+6,yy+2)
        p3 = (xx,yy-3)
        drawer.polygon([p1,p2,p3],outline=out,fill=col)
        p0 = (xx-3,yy-5)
        p1 = (xx-3,yy+3)
        p2 = (xx+3,yy+3)
        p3 = (xx+3,yy-5)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx-1,yy-9)
        p1 = (xx+1,yy-9)
        p2 = (xx+1,yy+4)
        p3 = (xx-1,yy+4)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
    def drawCapital(self,drawer,xx,yy,col,out):
        p1 = (xx-6,yy+2)
        p2 = (xx+6,yy+2)
        p3 = (xx,yy+9)
        drawer.polygon([p1,p2,p3],outline=out,fill=col)
    def drawSelf(self,drawer):
        col = self.culture.bannerColor
        out = (0,0,0)
        if self.culture.origin == self.node:
            self.drawCapital(drawer,self.node.x,self.node.y,col,out)
        if self.population <= self.popThresholds[0]:
            self.drawTent(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= self.popThresholds[1]:
            self.drawHut(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= self.popThresholds[2]:
            self.drawVillage(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= self.popThresholds[3]:
            self.drawTown(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= self.popThresholds[4]:
            self.drawCity(drawer,self.node.x,self.node.y,col,out)
        else:
            self.drawMetropolis(drawer,self.node.x,self.node.y,col,out)
    def popNotes(self):
        return self.cityNotes()
    def cityNotes(self):
        s = self.name + " ("
        s += self.cType(self.population) + ")\n\n"
        s += "Governance: " + self.governanceName() + "\n\n"
        s += "Civilian population: " + str(self.population) + "\n\n"
        s += "Number of notable entities at this location: " + str(len(self.node.entities)) + "\n\n"
        return s

class Road:
    def __init__(self,n0,n1,c):
        self.built = 0
        self.city = c
        self.city.roads.append(self)
        self.city.culture.myMap.roads.append(self)
        self.start = n0
        self.destination = n1
        self.current = n0
        self.nodes = []
        self.nodes.append(self.start)
        self.size = 2
        self.color = (20,120,104)
    def build(self):
        if self.built == 1:
            return -1
        currDist = 1000000000
        nextPick = None
        for n in self.current.neighbors:
            nDist = n.dist(self.destination)
            if nDist < currDist:
                if n.watery() == 0:
                    nextPick = n
                    currDist = nDist
        if nextPick == None:
            return -1
        self.linkRoad(self.current,nextPick)
        self.current = nextPick
        if self.current == self.destination:
            self.built = 1
            self.destination.city.roads.append(self)
    def linkRoad(self,n0,n1):
        self.nodes.append(n1)
        if n1 not in n0.roads:
            n0.roads.append(n1)
        if n0 not in n1.roads:
            n1.roads.append(n0)
    def drawSelf(self,drawer):
        nds = [x.coords() for x in self.nodes]
        drawer.line(nds,fill=self.color,width=self.size)

class ResourceRegion:
    def __init__(self,c,m):
        self.rootCity = c
        self.myMap = m
        self.myMap.resourceRegions.append(self)
        self.culture = self.rootCity.culture
        self.nodes = []
        self.addNode(self.rootCity.node)
        self.resources = [0,0]
        self.updateReg()
    def expungeReg(self):
        self.myMap.resourceRegions.remove(self)
    def addNode(self,node):
        if node.resourceRegion != None:
            node.resourceRegion.nodes.remove(node)
        node.resourceRegion = self
        self.nodes.append(node)
    def cityCount(self):
        c = 0
        q = 0
        for p in self.nodes:
            if p.city != None:
                c += 1
                q += p.city.population
        self.totalCities = c
        self.totalPop = q
    def sumResources(self):
        self.resources = [0,0] # [Food resources, industrial resources]
        rawPlant = 0
        rawMetal = 0
        rawAnimal = 0
        for p in self.nodes:
            if p.landmass != None:
                rawPlant += p.vegetation
                rawMetal += p.metallicity
                rawAnimal += p.herbivores + (p.carnivores/3)
            else:
                rawAnimal += 0.17
                rawPlant += 0.1
                rawMetal += 0.02
        m = self.culture.value.mainValues
        if "metallurgists" in m:
            rawMetal = rawMetal*1.25
        if "agriculture" in m:
            rawPlant *= 1.15
        if "warriors" in m:
            rawAnimal *= 1.2
        if "simplicity" in m:
            rawMetal *= 0.75
            rawAnimal *= 0.75
            rawPlant *= 0.75
        rawMetal *= self.culture.tech["metallurgy"]
        rawPlant *= self.culture.tech["agriculture"]
        rawAnimal*= self.culture.tech["weaponry"]
        scale = self.myMap.resourceScale
        self.resources[0] = scale*(rawPlant + (rawAnimal*0.7))
        self.resources[1] = scale*(rawMetal + (rawAnimal*0.3) + (rawPlant*0.4))
        if "craftsmen" in m:
            self.resources[1] *= 1.05
        if "builders" in m:
            self.resources[1] *= 1.05
        if "collectivists" in m:
            self.resources[0] *= 1.05
    def updateReg(self):
        for p in self.nodes:
            if p.city != None:
                p.resourceDist = math.log(p.city.population)
        for p in self.nodes:
            for k in p.neighbors:
                if k.resourceRegion != self and k.resourceRegion != None:
                    # In this case, it's a different region OTHER than noRegion.
                    if k.resourceRegion.culture != self.culture:
                        # Interact with a different region of a different society...
                        # This is a placeholder. In the future there will probably be wars here.
                        k.resourceRegion = k.resourceRegion
                    elif k.resourceRegion.culture == self.culture:
                        # Interacting with a different region of the same society...
                        if k.resourceDist < p.resourceDist:
                            self.addNode(k)
                            k.resourceDist = (p.resourceDist-1)/2
                elif k.resourceRegion == None:
                    if k.resourceDist < p.resourceDist:
                        self.addNode(k)
                        k.resourceDist = (p.resourceDist-1)/2
                elif k.resourceRegion == self:
                    if k.resourceDist < p.resourceDist:
                        k.resourceDist = (p.resourceDist-1)/2
        self.cityCount()
        self.sumResources()
        if len(self.nodes) == 0:
            self.expungeReg()
            return -1

class Culture:
    def __init__(self,n,m):
        self.origin = n
        self.myMap = m
        self.myMap.cultures.append(self)
        self.influences = Influence(self.myMap,self.origin,1)
        self.value = Values(self.myMap,self.influences)
        self.society = self.setSociety()
        self.setGenderBalance()
        self.language = Language(self)
        self.name = self.language.name
        self.populations = {}
        self.magic = []
        self.favoriteDeity = None
        self.generateMythology()
        self.title = self.setTitle()
        self.flag = Flag(self)
        self.bannerColor = self.flag.colors[0]
        self.leaderTitle = self.setLeaderTitle()
        self.cities = []
        if self.leaderCount == 1:
            self.ppp = "person"
        else:
            self.ppp = "group"
        self.leader = Population(self,t=self.leaderTitle,i=random.randint(5,25),p=self.leaderCount,kind=self.ppp,node=self.origin,prf="politician")
        self.totalPop = self.populationCount()
        self.oldAge = 65
        self.tech = {}
        self.items = []
        for k in list(self.myMap.technologies.keys()):
            self.tech[k] = 1+(0.1*random.random())
        for p in range(50):
            self.updateTech()
        self.unitcount = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        self.cultureFace = None
        self.cultureOpinions = {}
    def happiness(self):
        p = 0
        h = 0
        for f in self.cultureCities():
            p += f.population
            h += f.population*f.happiness()
        h = h/p
        return h
    def populationCount(self):
        p = 0
        for f in self.cultureCities():
            p += f.population
        return p
    def closestCities(self,otherCulture):
        ours = self.cultureCities()
        theirs = otherCulture.cultureCities()
        ourCity = random.choice(ours)
        theirCity = random.choice(theirs)
        dist = ourCity.distanceToCity(theirCity)
        for c in ours:
            for g in theirs:
                d = c.distanceToCity(g)
                if d < dist:
                    dist = d
                    ourCity = c
                    theirCity = g
        return [ourCity,theirCity]
    def cultureCities(self):
        cc = []
        for f in self.myMap.cities:
            if f.culture == self:
                cc.append(f)
        return cc
    def initOpinions(self):
        for c in self.myMap.cultures:
            if not (c.name == self.name):
                newOpinion = Opinion(self,c)
                self.cultureOpinions[c.name] = newOpinion
    def updateOpinions(self):
        for o in self.cultureOpinions:
            self.cultureOpinions[o].updateOpinion()
    def updateTech(self):
        for t in self.tech.keys():
            multiplier = 0.005 + random.choice([-0.001,0.000,0.001])
            m = self.value.mainValues
            if "metallurgists" in m:
                if t == "metallurgy":
                    multiplier = multiplier*1.2
            if "builders" in m:
                if t == "defense":
                    multiplier = multiplier*1.2
            if "collectivists" in m:
                if t == "government":
                    multiplier = multiplier*1.4
            if "agriculture" in m :
                if t == "agriculture":
                    multiplier = multiplier*1.4
            if "materialists" in m :
                if t == "research":
                    multiplier = multiplier*1.2
            if "warriors" in m :
                if t == "weaponry":
                    multiplier = multiplier*1.2
            if "craftsmen" in m:
                if t == "production":
                    multiplier = multiplier*1.4
            if "freedom" in m:
                if t == "equality":
                    multiplier = multiplier*1.2
            if "worshippers" in m or "astrology" in m:
                if t == "art" or t == "philosophy":
                    multiplier = multiplier*1.3
            if "simplicity" in m:
                multiplier = multiplier*0.75
            if "tribe" in self.society.lower():
                multiplier = multiplier*0.4
            if "shaman" in self.society.lower():
                multiplier = multiplier*0.75
            if "hunter" in self.society.lower():
                multiplier = multiplier*0.4
            if self.society == "Scholars":
                multiplier = multiplier*1.1
                if t == "research" or t == "philosophy":
                    multiplier = multiplier*1.4
            if self.society == "Blacksmiths":
                if t == "metallurgy":
                    multiplier = multiplier*1.2
            if "artisan" in self.society:
                if t == "production":
                    multiplier = multiplier*1.2
            if "agricultur" in self.society.lower():
                if t == "agriculture":
                    multiplier = multiplier*1.3
            if t == "research":
                multiplier = multiplier*(self.tech["philosophy"]**0.2)
            multiplier = multiplier*(self.tech["research"]**0.3333)
            multiplier = multiplier-(self.tech[t]*0.0005)
            self.tech[t] = self.tech[t]*(multiplier+1)
    def updatePops(self):
        for f in list(self.populations.keys()):
            if f in self.populations:
                self.populations[f].agePop(self.myMap.timeScale)
        if self.leader == None or self.leader.number == 0:
            pp = self.leader
            parens = []
            if self.society in ["Hegemony","Empire","Imperium","Monarchy"]:
                for offspring in pp.kids:
                    if offspring.dead == 0 and self.leader == pp:
                        self.leader = offspring
                        self.leader.title = self.leaderTitle
                        self.leader.node = self.origin
            if self.leader == pp:
                self.leader = Population(self,t=self.leaderTitle,i=random.randint(5,25),p=self.leaderCount,kind=self.ppp,node=self.origin,prf="politician",pars=parens)
        if self.leader.number < self.leaderCount:
            replaceAmount = self.leaderCount-self.leader.number
            currentContribution = math.floor(self.leader.age*(self.leader.number/self.leaderCount))
            newContribution = math.floor((random.uniform(26,60)*(replaceAmount/self.leaderCount)))
            self.leader.age = currentContribution+newContribution
            self.leader.number = self.leaderCount
        nn = 0.2
        pp = self.populationCount()
        figs = len(self.populations)
        pmod = clamp((pp**nn)/(0.90),0,10)*0.01
        redc = clamp(figs/75,0,1)*0.1
        while random.random() > ((0.985-pmod)+redc):
            ctnodes = [c.node for c in self.cultureCities()]
            fig = Population(self,p=1,kind="person",node=random.choice(ctnodes))
    def updateCulture(self):
        self.updateOpinions()
        self.updateTech()
        if self.deities[0].age % self.electionYear == 0:
            self.election()
    def election(self):
        if self.leaderCount > 1:
            self.leader.kill((1-self.happiness())*self.leader.number)
        else:
            if self.happiness() < 0.5 or random.random() > 0.65:
                ll = self.leader
                ll.title = ""
                if random.random() > (0.1*len(self.listOfProfessions("politician"))):  # Either elect a new politician
                    self.leader = Population(self,t=self.leaderTitle,p=self.leaderCount,kind=self.ppp,node=self.origin,prf="politician")
                else:   # Or pick from politicians they already have
                    politicians = self.listOfProfessions("politician")
                    newleader = random.choice(politicians)
                    newleader.travel(self.origin)
                    self.leader = newleader
        self.leader.title = self.leaderTitle
        ev = Event(self.myMap,a=1,kind="election",sub=self.leader,loc=self.origin)
    def listOfProfessions(self,kind,alive=1):
        picks = []
        for p in self.populations.keys():
            if self.populations[p].profession == kind:
                if (self.populations[p].dead == 0 or alive == 0):
                    picks.append(self.populations[p])
        return picks
    def setSociety(self):
        m = self.value.mainValues
        if "warriors" in m and "collectivists" in m and "worshippers" in m:
            return "Hegemony"
        if ("collectivists" in m and "worshippers" in m and 
            ("freedom" not in m and "materialists" not in m)):
            return "Monarchy"
        if "greed" in m and "worshippers" in m and "builders" in m:
            return "Empire"
        if "warriors" in m and "builders" in m and "worshippers" in m and "freedom" not in m and "traders" not in m:
            return "Nation-state"
        if ("warriors" in m and "greed" in m and "builders" in m and 
            ("freedom" in m or "travelers" in m or "sailors" in m)):
            return "Imperium"
        if ("materialists" in m and "astrology" in m and 
            ("superstition" not in m or "worshippers" not in m)):
            return "Astronomers"
        if ("freedom" in m and "materialists" in m and ("superstition" not in m and "shamans" not in m)):
            return "Scholars"
        if "freedom" in m and "collectivists" in m and "simplicity" in m:
            return "Paleolithic tribe"
        if "collectivists" in m and "agriculture" in m and "materialists" in m:
            return "Agricultural commune"
        if "collectivists" in m and "agriculture" in m and "simplicity" in m:
            return "Farming commune"
        if "worshippers" in m and "warriors" in m and ("superstition" in m or "collectivists" in m):
            return "Religious zealots"
        if "builders" in m and "collectivists" in m and "materialists" in m:
            return "Socialists"
        if "metallurgists" in m and "builders" in m and "craftsmen" in m and "traders" in m:
            return "Merchant artisans"
        if "freedom" in m and "greed" in m and "traders" in m:
            return "Liberal capitalists"
        if "builders" in m and "agriculture" in m and ("travelers" in m or "sailors" in m):
            return "Colonists"
        if ("craftsmen" in m or "metallurgy" in m or "builders" in m) and "simplicity" in m and ("naturalists" in m or "shamans" in m):
            return "Naturalist artisans"
        if ("craftsmen" in m or "metallurgy" in m or "builders" in m) and "simplicity" in m and ("superstition" in m or "astrology" in m):
            return "Traditionalist artisans"
        if (("travelers" in m or "sailors" in m) and "greed" in m and "warriors" in m):
            return "Raiders"
        if (("travelers" in m or "sailors" in m) and "greed" in m and "simplicity" in m):
            return "Scavengers"
        if "shamans" in m and "warriors" in m and ("astrology" in m or "superstition" in m or "worshippers" in m):
            return "Shamanistic warriors"
        if "travelers" in m and "simplicity" in m and "freedom" in m:
            return "Hunter-gatherer tribe"
        if "astrology" in m and "superstition" in m and "worshippers" in m:
            return "Religious sovereignty"
        if "warriors" in m and "collectivists" in m and ("travelers" in m or "sailors" in m):
            return "Revolutionary commune"
        if ("builders" in m and "metallurgists" in m and "craftsmen" in m and 
            "agriculture" in m and "collectivists" in m):
            return "Syndicalists"
        if ("shamans" in m and ("naturalists" in m or "astrology" in m or "simplicity") and 
            "superstition" in m):
            return "Shamanic tribe"
        if "metallurgists" in m and "builders" in m and "craftsmen" in m and "materialists" in m:
            return "Cooperative artisans"
        if "greed" in m and "sailors" in m and "warriors" in m:
            return "Pirates"
        if "freedom" in m and "traders" in m and ("builders" in m or "craftsmen" in m or "metallurgists" in m):
            return "Liberal merchant-artisans"
        if "collectivists" in m and "freedom" in m:
            return "Social democracy"
        if "builders" in m and "agriculture" in m and "traders" in m:
            return "Mercantile folk"
        if "materialists" in m and "metallurgists" in m and ("craftsmen" in m or "builders" in m):
            return "Blacksmiths"
        if "collectivists" in m and "worshippers" in m:
            return "Religious collective"
        if "traders" in m and "collectivists" in m:
            return "Co-operative"
        if "collectivists" in m and "agriculture" in m and "naturalists" in m:
            return "Agricultural naturalists"
        if "travelers" in m and "simplicity" in m and ("metallurgists" in m or "craftsmen" in m or "builders" in m):
            return "Nomadic artisans"
        if "agriculture" in m and "worshippers" in m:
            return "Religious agriculturalists"
        if ("travelers" in m or "sailors" in m) and ("traders" in m or "greed" in m) and "freedom" in m:
            return "Independent merchants"
        if (("astrology" in m and "superstition" in m) or 
            ("superstition" in m and "worshippers" in m) or
            ("woshippers" in m and "astrology" in m)):
            return "Religious sovereignty"
        if "travelers" in m and "simplicity" and ("naturalists" in m or "superstition" in m or "astrology" in m or "shamans" in m):
            return "Nomadic tribe"
        if "travelers" in m and "sailors" in m and "traders" in m:
            return "Mariners"
        if "traders" in m and "freedom" in m:
            return "Merchants"
        if "shamans" in m:
            return "Shamans"
        if "worshippers" in m and "superstition" in m:
            return "Religious sovereignty"
        if "traders" in m and "greed" in m:
            return "Capitalists"
        if "collectivists" in m:
            return "Communalists"
        if "simplicity" in m:
            return "Tribe"
        if "freedom" in m:
            return "Liberals"
        if "metallurgists" in m:
            return "Blacksmiths"
        if "agriculture" in m:
            return "Agriculturalists"
        if "craftsmen" in m:
            return "Craftsmen"
        if "travelers" in m:
            return "Nomads"
        if "greed" in m and "warriors" in m:
            return "Raiders"
        if "sailors" in m and "warriors" in m:
            return "Pirates"
        if "traders" in m:
            return "Traders"
        if "greed" in m:
            return "Capitalists"
        return "Mixed society"
    def setTitle(self):
        titleDict = dict.fromkeys(["Nation-state"],["Nation"])
        titleDict.update(dict.fromkeys(["Religious sovereignty",
                      "Religious zealots",
                      "Religious agriculturalists",
                      "Religious collective"],["Theocracy","Ecclesiarchy","Order","Caliphate","Cult","See","Papacy"]))
        titleDict.update(dict.fromkeys(["Agriculturalists",
                       "Farming commune",
                       "Agricultural commune",
                       "Agricultural naturalists"],["Farmers","Yeomen","Peasants","Cultivators"]))
        titleDict.update(dict.fromkeys(["Empire",
                        "Hegemony",
                        "Imperium"],["Empire"]))
        titleDict.update(dict.fromkeys(["Monarchy"],["Kingdom"]))
        titleDict.update(dict.fromkeys(["Nomadic artisans",
                        "Nomads",
                        "Scavengers",
                        "Nomadic tribe"],["Nomads","Pilgrims","Travelers","Bands"]))
        titleDict.update(dict.fromkeys(["Liberal capitalists",
                         "Liberal merchant-artisans",
                         "Merchant artisans",
                         "Traders",
                         "Independent merchants",
                         "Mercantile folk",
                         "Merchants",
                         "Capitalists",
                         "Colonists",
                         "Mariners"],["Caravans","Proprietors","Holdings","Shipping","Commerce","Traders","Merchants"]))
        titleDict.update(dict.fromkeys(["Blacksmiths",
                          "Traditionalist artisans",
                          "Naturalist artisans",
                          "Cooperative artisans",
                          "Craftsmen"],["Artisans","Craftsmen","Guild"]))
        titleDict.update(dict.fromkeys(["Socialists",
                           "Syndicalists",
                           "Revolutionary commune",
                           "Communalists",
                           "Co-operative"],["People's Union","Union","Collective"]))
        titleDict.update(dict.fromkeys(["Shamanistic warriors",
                            "Shamanic tribe",
                            "Shamans"],["Mystics","Shamanate","Order"]))
        titleDict.update(dict.fromkeys(["Pirates",
                             "Raiders"],["Brigands","Raiders","Legion","Bands"]))
        titleDict.update(dict.fromkeys(["Social democracy",
                              "Liberals"],["Republic"]))
        titleDict.update(dict.fromkeys(["Scholars",
                                        "Astronomers"],["Institute","Academy","College","Order"]))
        t = "People"
        for h in titleDict.keys():
            if h == self.society:
                t = random.choice(titleDict[h])
        self.electionYear = 1000000000
        if self.society in ["Socialists","Syndicalists","Revolutionary commune","Communalists","Co-operative"]:
            self.electionYear = random.randint(1,5)
        if self.society in ["Social democracy","Liberals","Scholars","Astronomers"]:
            self.electionYear = random.randint(3,10)
        if "commune" in self.society:
            self.electionYear = random.randint(1,5)
        if "liberal" in self.society:
            self.electionYear = random.randint(3,10)
        if "cooperative" in self.society:
            self.electionYear = random.randint(8,20)
        return t
    def setLeaderTitle(self):
        self.leaderCount = 1
        s = ""
        s2 = ""
        if self.society == "Nation-state":
            return (random.choice(["Supreme ","High ","Lord ",""]) +
                    random.choice(["Commissioner","Chancellor","Harbinger"]))
        if (self.society == "Religious sovereignty" or self.society == "Religious zealots" or
            self.society == "Religious agriculturalists"):
            return (random.choice(["Grand ","High ","Supreme ","Holy "]) +
                    random.choice(["Pontiff","Priest","Shepherd"]))
        if (self.society == "Agriculturalists" or self.society == "Farming commune" 
            or self.society == "Agricultural commune" or self.society == "Agricultural naturalists"):
            self.leaderCount = random.choice([1,random.randint(2,10)])
            if self.leaderCount == 1:
                s = ""
            else:
                s = random.choice(["Council","Assembly","Soviet","Conference","Directorate"]) + " of "
                s2 = "s"
            return (s + random.choice(["Head ","Chief ","Master "])+
                    random.choice(["Farmer","Agronomist","Foreman"])+s2)
        if self.society == "Imperium" or self.society == "Hegemony" or self.society == "Empire":
            return "Emperor"
        if self.society == "Monarchy":
            return "Monarch"
        if self.society == "Nomadic artisans" or self.society == "Nomads" or self.society == "Scavengers":
            return (random.choice(["Chief ","Head ","Elder "]) +
                    random.choice(["Captain","Dignitary","Herald"]))
        if (self.society == "Liberal capitalists" or self.society == "Liberal merchant-artisans" or self.society == "Merchant artisans" or 
            self.society == "Traders" or self.society == "Independent merchants" or self.society == "Mercantile folk"
            or self.society == "Merchants" or self.society == "Capitalists" or self.society == "Colonists" or self.society == "Mariners"):
            self.leaderCount = random.choice([1,random.randint(2,10)])
            if self.leaderCount > 1:
                s = random.choice(["Cabinet","Assembly","Board","Committee"]) + " of "
                s2 = "s"
            return (s + random.choice(["Primary ","Head ","Chief ",""]) +
                    random.choice(["Executive","Director","Superintendent"]) + s2)
        if (self.society in ["Blacksmiths","Traditionalist artisans","Naturalist artisans","Cooperative artisans","Craftsmen"]):
            self.leaderCount = random.choice([1,random.randint(2,10)])
            if self.leaderCount > 1:
                s = random.choice(["Council","Assembly","Congress"]) + " of "
                s2 = "s"
            return (s + random.choice(["Master ","Elder ","Grandmaster "]) +
                    random.choice(["Artificer","Builder","Craftsperson"]) + s2)
        if (self.society == "Socialists" or self.society == "Syndicalists" or self.society == "Revolutionary commune"
            or self.society == "Communalists" or self.society == "Co-operative"):
            self.leaderCount = random.choice([1,random.randint(2,20)])
            if self.leaderCount > 1:
                s = random.choice(["Council","Assembly","Soviet","Conference","Directorate"]) + " of "
                s2 = "s"
            return (s+random.choice(["Prime ","Chief ","Central ",""]) +
                    random.choice(["Director","Governer","Speaker","Chairperson"])+s2)
        if (self.society == "Shamanistic warriors" or self.society == "Shamanic tribe"
            or self.society == "Shamans"):
            return (random.choice(["Elder","High","Grand","Ancestral"])+ " " +
                    random.choice(["Medicine Man","Seer","Shaman"]))
        if self.society == "Pirates" or self.society == "Raiders":
            return (random.choice(["Chief ","Head ",""]) +
                    random.choice(["Captain","Commander","Warlord"]))
        if self.society == "Scholars" or self.society == "Astronomers":
            self.leaderCount = random.choice([1,random.randint(2,10)])
            if self.leaderCount > 1:
                s = random.choice(["Council","Assembly","Congress"]) + " of "
                s2 = "s"
            return (s + random.choice(["Master ","Elder ","Grandmaster ",""]) +
                    random.choice(["Dean","Chancellor","Professor"]) + s2)
        if self.society == "Social democracy" or self.society == "Liberals":
            self.leaderCount = random.choice([1,1,random.randint(2,10)])
            if self.leaderCount > 1:
                s = random.choice(["Congress","Chamber","Parliament","Ministry","Senate"]) + " of "
                s2 = "s"
            return s + random.choice(["President","Speaker","Minister","Representative","Premier"]) + s2
        return "Chief"
    def setGenderBalance(self):
        man = math.floor(random.uniform(0,100))
        woman = math.floor(random.uniform(0,100))
        nb = math.floor(random.uniform(0,100))
        self.genderSpread = [1 for n in range(man)]+[2 for n in range(woman)]+[3 for n in range(nb)]
    def generateMythology(self):
        self.deities = []
        m = self.value.mainValues
        tiers = 4
        maxtiers = 6
        if "superstition" in m:
            t = random.randint(1,3)
            tiers += t
        if "astrology" in m:
            t = random.randint(0,1)
            tiers += t
        if "worshippers" in m:
            t = random.randint(1,2)
            tiers += t
        if "shamans" in m:
            t = random.randint(1,2)
            tiers += t
        if "naturalists" in m:
            t = random.randint(0,1)
            tiers += t
        if "materialists" in m:
            t = -1
            tiers += t
        if "simplicity" in m:
            t = -1
            tiers += t
        if (("superstition" not in m) and ("worshippers" not in m) and 
            ("astrology" not in m) and ("shamans" not in m) and ("naturalists" not in m)):
            t = -1
            tiers += t
        tiers = clamp(tiers,0,maxtiers)
        titles = ["","","","","","","",""]
        k = 0
        while k < maxtiers:
            title = ""
            if k == 0:
                title += random.choice(["foundational ","formless ","empty ",
                                        "primordial ","chaotic ","cosmic ",
                                        "astral ","metaphysical ","infinite ",""])
            if k == 1:
                title += random.choice(["primordial ","primal ","primeval ","infinite ",
                                        "omnipotent ","universal ","transcendental ",""])
            if k == 2:
                title += random.choice(["elder ","old ","ancient ","celestial ",
                                        "cosmic ","nameless ","unspeakable ",""])
            if k == 3:
                title += random.choice(["great ","vast ","titanic ","exalted ",
                                        "divine ","ineffable ",""])
            if k == 4:
                title += random.choice(["major ","immortal ","astral ",
                                        "empyrean ","empyreal "])
            if k == 5:
                title += random.choice(["minor ","mortal ","ascetic ","young ","living ",""])
            title = title
            if k == 0:
                title += random.choice(["universe","void","expanse","firmament","reality"])
            else:
                title += random.choice(["titan","god","being","leviathan",
                                        "entity","deity","godhead","overbeing"])
            titles[k] = title
            if k == 0:
                k = clamp(maxtiers-tiers,1,maxtiers)
            else:
                k += 1
        aa = random.random()*8
        ageOf = round(self.myMap.age*aa)
        self.mythAge = ageOf+5
        rr = 1.2
        f = Population(self,n=self.language.genName(),t=titles[0],a=ageOf,kind="location",i=random.randint(30,110))
        f.gender = 0
        if ageOf < self.myMap.age:
            e = Event(m=self.myMap,a=ageOf,kind="genesis",sub=f)
            e.importance = 100
        self.deities.append(f)
        for k in range(clamp(maxtiers-tiers,1,maxtiers),tiers):
            entities = math.ceil((random.randint(0,2)+k)/2)
            if random.random() < 0.3:
                entities = 0
            if random.random() < 0.25:
                entities = random.randint(5,7)
            for e in range(entities):
                mm = 1
                if len(self.deities) >= 2:
                    mm = 2
                if random.random() < (0.43 + (k*0.1)):
                    parent = random.randint(1,mm)
                else:
                    parent = 0
                age = round(ageOf/(len(self.deities)*rr))
                if age > self.myMap.age:
                    parent = 0
                nom = self.language.genName()
                if random.random() < (0.25 + (k*0.1)):
                    nom += " " + self.language.genName()
                ei = math.floor(90/k)+random.randint(5,20)
                ent = Population(self,n=nom,t=titles[k],a=age,kind="deity",i=ei,pars=[])
                ent.kids = []
                for u in range(parent):
                    pp = random.choice(self.deities)
                    while pp in ent.parents:
                        pp = random.choice(self.deities)
                    ent.parents.append(pp)
                    pp.kids.append(ent)
                ent.associate()
                if random.random() < 0.2/k:
                    newSpell = Magic(ent)
                    newSpell.strength = clamp(1/k,0.5,1)
                self.deities.append(ent)
                if age < self.myMap.age:
                    e = Event(m=self.myMap,a=age,kind="birth",sub=ent,actrs=ent.parents)
                    e.importance = ent.importance*random.uniform(0.66666,1.33333)
                    ent.birthEvent = e
        self.favoriteDeity = random.choice(self.deities)
    def generateCultureFace(self,mode):
        filename = "./generated/face_" + self.name + ".gif"
        if self.cultureFace == None:
            res = 192
            self.cultureFace = Face(self,x=res)
            self.cultureFace.generateCultureFace()
            img = Image.new('HSV',(res,res),(255,0,255))
            drawer = ImageDraw.Draw(img)
            self.cultureFace.drawSelf(drawer)
            img = img.convert('RGB')
            img.save(filename,"GIF")
        if mode == 0:
            return filename
        else:
            return self.cultureFace
    def nameOfCapital(self):
        capitalName = "";
        synonymSeed = seedNum(self.name)
        if self.society in ["Monarchy","Imperium","Empire","Hegemony"]:
            capitalName += synonym("palace",synonymSeed,0).capitalize()
            if synonymSeed % 2 == 0:
                capitalName += " of " + self.leader.name[1].title()
            else:
                capitalName += " of " + self.name.title()
        #
        #
        elif self.society in ["Religious sovereignty",
                      "Religious zealots",
                      "Religious agriculturalists",
                      "Religious collective"] or "shaman" in self.society.lower():
            capitalName += synonym("cathedral",synonymSeed,0).title()
            if synonymSeed % 2 == 0:
                capitalName += " of " + self.favoriteDeity.justName()
            else:
                capitalName += " of " + self.name.title() + " "
        #
        #
        elif "tribe" in self.society.lower():
            capitalName += synonym("longhouse",synonymSeed,0).title()
            if synonymSeed % 3 == 0:
                capitalName += " of " + self.leader.name[1].title()
            elif synonymSeed % 3 == 1:
                capitalName += " of " + self.name
        #
        #
        elif self.leaderCount > 1:
            if synonymSeed % 2 == 0:
                capitalName += "National "
            else:
                capitalName += self.name.capitalize() + " "
            capitalName += synonym("parliament",synonymSeed,0).title()
        #
        #
        else:
            if synonymSeed % 2 == 0:
                capitalName += "National "
            capitalName += synonym("office",synonymSeed,0).title()
            if synonymSeed % 2 == 1:
                capitalName += " of " + self.name.title()
        return capitalName
    def nameOfTownHall(self,townName):
        townName = townName.title()
        capitalName = "";
        synonymSeed = seedNum(self.name)
        if self.society in ["Monarchy","Imperium","Empire","Hegemony"]:
            capitalName += synonym("office",synonymSeed,0).title()
            if synonymSeed % 2 == 0:
                capitalName += " of " + self.leader.name[1].title()
            else:
                capitalName += " of " + self.name.title()
        #
        #
        elif self.society in ["Religious sovereignty",
                      "Religious zealots",
                      "Religious agriculturalists",
                      "Religious collective"] or "shaman" in self.society.lower():
            capitalName += synonym("church",synonymSeed,0).title()
            if synonymSeed % 2 == 0:
                capitalName += " of " + self.favoriteDeity.justName()
            else:
                capitalName += " " + townName
        #
        #
        elif "tribe" in self.society.lower():
            capitalName += synonym("longhouse",synonymSeed,0).title()
            if synonymSeed % 3 == 0:
                capitalName += " of " + townName
            elif synonymSeed % 3 == 1:
                capitalName += " of " + self.name
        #
        #
        elif self.leaderCount > 1:
            if synonymSeed % 2 == 0:
                capitalName += "Local "
            else:
                capitalName += townName + " "
            capitalName += synonym("governance",synonymSeed,0).title()
        #
        #
        else:
            if synonymSeed % 2 == 0:
                capitalName += "National "
            else:
                capitalName += self.name.capitalize() + " "
            capitalName += synonym("office",synonymSeed,0).title()
        return capitalName
    def shortName(self):
        name = ""
        name += self.name + " " + self.title
        return name
    def information(self):
        info = ""
        info += self.name +" "+ self.title + "\n"
        info += "("+self.society+")" + "\n"
        return info
    def mythNotes(self):
        s = ""
        for d in self.deities:
            if self.deities.index(d) != 0:
                s += "\n\n"
            if d.age > self.myMap.age:
                s += "Since the beginning of time, "
            else:
                s += str(nearestHundred(d.age)) + " years ago, "
            s += "the " + d.nameFull()
            if d.age > self.myMap.age:
                s += " has always existed"
            else:
                s += " was born from "
                if len(d.parents) == 0:
                    s += "nothing"
                elif len(d.parents) == 1:
                    q = d.parents[0]
                    s += "the " + q.nameFull()
                else:
                    q0 = d.parents[0]
                    q1 = d.parents[1]
                    s += "the " + q0.nameFull()
                    s += " and the " + q1.nameFull()
            s += "."
        return s
    def cultureNotes(self):
        s = self.name + " " + self.title + "\n"
        s += "Society type: " + self.society + "\n\n"
        if self.leaderCount == 1:
            s += "Leader: The " + self.leader.nameFull() + "\n\n"
        else:
            s += "Leading body: The " + self.leader.nameFull() + "\n\n"
        s += "Civilian population: " + str(self.populationCount()) + "\n\n"
        s += "Capital: " + self.origin.city.name + "\n"
        return s
    def techNotes(self):
        s = ""
        for k in self.tech.keys():
            lv = self.tech[k]
            lvl = techTier(lv)
            testlvl = abs(round(lvl))
            lvldiff = abs(testlvl-lvl)
            if lvl > testlvl:
                if lvldiff > 0.32:
                    ss = " > "
                elif lvldiff > 0.13:
                    ss = " >= "
                else:
                    ss = " ~= "
            else:
                if lvldiff > 0.32:
                    ss = " < "
                elif lvldiff > 0.13:
                    ss = " <= "
                else:
                    ss = " ~= "
            lvl = math.floor(lvl)
            s += "" + k
            s += "" + ss
            s += self.myMap.techtiers[lvl] + " humanity. ("
            s += str(lv) + ")"
            s += "\n\n"
        return s
    def drawPops(self,drawer):
        for p in self.populations.keys():
            self.populations[p].drawSelf(drawer)

class Opinion:
    def __init__(self,c,o):
        self.culture = c
        self.other = o
        # 0 - no knowledge; 1 - only know of each other; 2 - contact.
        self.knowledge = 0
        # [1,0,0] Friendly
        # [-1,0,0] Hostile
        # [0,1,0] Stronger
        # [0,-1,0] Weaker
        # [0,0,1] Interventionist
        # [0,0,-1] Isolationist
        self.status = [0,0,0]
        self.activeActions = []
    def updateOpinion(self):
        actions = {}
        # last item in the list is the radius of the action
        actions["war"] = [-1,1,1,0.75]
        actions["alliance"] = [1,-0.5,1,0.7]
        actions["aid"] = [1,1,1,0.8]
        actions["close borders"] = [-0.85,-0.5,-1,0.75]
        actions["non aggression pact"] = [1,-0.5,-1,1]
        actions["sanction"] = [-0.5,0.5,-1,0.75]
        if self.knowledge == 0:
            return
        return 1

# This can represent either an entity, a group, OR sometimes an imaginary location (in the case of mythology)
class Population:
    #                 culture, name, title, age, count, kind, node, profession, importance, parents
    def __init__(self,c=None,n=None,t="",a=None,p=1,kind="person",node=None,prf=None,i=None,pars=[]):
        self.tt = "pop"
        self.dead = 0
        self.parents = pars
        self.kids = []
        self.culture = c
        self.number = p
        self.profession = prf
        self.kind = kind
        self.speed = 1.0
        self.importance = math.floor((random.random()**2)*23)
        self.magic = []
        self.birthEvent = None
        self.deathEvent = None
        if i != None:
            self.importance = i
        if len(self.parents) > 0:
            lst = [p.importance for p in self.parents]
            self.importance = self.importance + (0.2*(sum(lst)/len(lst)))
        self.baseImportance = self.importance
        self.measure = "tall"
        self.measurement = random.uniform(1.45,2)
        if (prf == None and (self.kind == "group" or self.kind == "person")):
            self.professions = ["artist","artist","artist",
                                "philosopher","researcher",
                           "politician","engineer","tactician",
                           "doctor","farmer","sociologist",
                           "blacksmith","blacksmith","blacksmith",
                           "artist","artist","philosopher",
                           "researcher","tactician","explorer","explorer","magician","magician","magician"]
            # Here is where we can add/remove items from the list to make certain professions
            # more or less common in certain societies.
            self.profession = random.choice(self.professions)
        self.field = None
        if self.profession != None:
            self.fields = {"artist":"art","philosopher":"philosophy",
                      "researcher":"research","politician":"government",
                      "engineer":"production","tactician":"weaponry",
                      "doctor":"medicine","farmer":"agriculture",
                      "sociologist":"equality","blacksmith":"metallurgy",
                      "craftsman":"production","roadbuilder":"transportation",
                      "explorer":"transportation","magician":"magic"}
            if self.kind != "army":
                self.field = self.fields[self.profession]
        else:
            self.field = None
        if n == None or n in self.culture.populations.keys():
            self.name = (self.culture.language.genName(),self.culture.language.genName())
            if self.number > 1:
                self.name = (self.culture.language.genName(),"")
        else:
            self.name = (n,"")
        if a == None:
            self.age = random.randint(23,50)
        else:
            self.age = a
        if t != "":
            self.title = t + " "
        else:
            self.title = t
        self.fullName = self.nameFull()
        self.culture.populations[self.name] = self
        self.pronouns = ["it","he","she","they"]
        self.possessive = ["its","his","her","their"]
        self.toBe = ["is","is","is","are"]
        self.gender = random.choice([1,2])
        if self.kind == "person":
            self.gender = random.choice(self.culture.genderSpread)
        if self.kind == "group" or self.kind == "army":
            self.gender = 0
        self.description = ""
        self.associations = []
        self.location = None
        if node != None:
            self.travel(node)
        self.terrain = 2    # 0=land, 1=water, 2=air(both)
        if self.kind == "beast":
            self.genBeast()
        self.immortal = 0
        if (self.kind == "deity" or self.kind == "beast"
            or self.kind == "location" or self.kind == "army" or self.kind == "group"):
            self.immortal = 1
        self.dead = 0
        if self.kind not in ["army","deity","beast","location","group"]:
            e = Event(m=self.culture.myMap,a=self.age,kind="birth",sub=self,actrs=self.parents,loc=self.location)
            e.importance = self.importance*random.uniform(0.6,2)
            self.birthEvent = e
        self.works = []
        self.units = ["assault infantry","guard infantry","siege","artillery","ranged infantry","mechanized","cavalry","beast"]
        self.rps = {"assault infantry":["artillery","siege"],
                    "guard infantry":["assault infantry","ranged infantry","cavalry"],
                    "siege":["guard infantry","artillery"],
                    "artillery":["guard infantry","assault infantry","ranged infantry","artillery","mechanized"],
                    "ranged infantry":["assault infantry","cavalry"],
                    "mechanized":["siege","cavalry","artillery","assault infantry"],
                    "cavalry":["assault infantry","artillery","siege","cavalry"]}
        self.unitbalance = {"assault infantry":[1,1,1],
                            "guard infantry":[1,1.5,0.8],
                            "siege":[1.25,0.75,2],
                            "artillery":[2,0.5,1.5],
                            "ranged infantry":[1.25,0.75,1.1],
                            "mechanized":[1.5,1.5,2],
                            "cavalry":[1.25,1,1.5]}
        if self.kind == "army":
            self.power = [1,1]
            self.gender = 0
            if self.profession == None:
                self.profession = "assault infantry"
            uindex = self.units.index(self.profession)
            self.culture.unitcount[uindex] += 1
            unitnumber = self.culture.unitcount[uindex]
            ss = self.culture.name.capitalize() + " " + ordinal(unitnumber) + " "
            q = 0
            if (self.profession in ["assault infantry","ranged infantry","siege"]):
                q = 1
            ss += synonym(self.profession,seed=seedNum(self.culture.name),exclusive=q)
            ss = string.capwords(ss)
            self.name = (ss,"")
        self.face = None
        self.path = None
    def associate(self):
        k = random.choice(list(self.culture.myMap.spheres))
        self.associations.append(k)
        self.descrip()
    def genBeast(self):
        self.gender = random.choice([1,2])
        humanoids = ["giant","elemental","demon","angel","hobgoblin"
                    "banshee","wendigo","naga","centaur","golem","werewolf",
                    "dryad","minotaur","ogre","skeleton","vampire",
                    "cyclops","gorgon","titan","troll","treefolk"]
        bipedalA = ["tyrannosaur","ostrich","emu"]
        bipedalB = ["chimpanzee","monkey","sloth","kangaroo","gorilla"]
        quadrupedsA = ["tiger","lion","lizard","wolf","rhinoceros",
                       "hyena","dingo","alligator","bear","crocodile",
                       "hydra","hippopotamus","basilisk"]
        quadrupedsB = ["chameleon","cattle","zebra","horse","elephant","armadillo",
                       "rabbit","antelope","deer","rat","mastodon","platypus","giraffe",
                       "tortoise","otter","frog","toad","brotosaur","boar",
                       "unicorn","goat","sheep","elk","anteater"]
        unconventional = ["snake","slug","roach","ant","spider","snail",
                          "slime","worm","scorpion","wurm","centipede","millipede"]
        sky = ["dragon","vulture","seagull","hawk","eagle","manticore","swan",
                  "duck","wasp","bumblebee","raven","bird","dragonfly","bat",
                  "hippogriff","pegasus","imp","albatross","sphinx"]
        ocean = ["whale","shark","jellyfish","squid","octopus","bass","trout","carp",
                   "koi","snake","snail","slug","seal","turtle","porpoise","penguin",
                   "otter","manatee","salamander","angler","dolphin","crab","cuttlefish",
                   "skate","ray","serpent","seahorse","nautilus","kraken","eel","angler",
                   "lobster"]
        if self.location.watery() == 1:
            self.terrain = 1
        if random.random() < 0.25:
            self.terrain = 2
        self.aggression = random.choice(["docile","neutral","aggressive"])
        if self.terrain == 0:
            if self.aggression == "docile":
                spheres = bipedalB+quadrupedsB+unconventional
            elif self.aggression == "neutral":
                spheres = bipedalB+quadrupedsB+unconventional+humanoids
            else:
                spheres = bipedalA+bipedalB+quadrupedsA+quadrupedsB+unconventional+humanoids
        elif self.terrain == 1:
            spheres = ocean
        else:
            spheres = sky
        self.species = random.choice(spheres)
        self.size = random.choice(["large","huge","gigantic"])
        if self.size == "large":
            bigness = 1
        elif self.size == "huge":
            bigness = 2
        elif self.size == "gigantic":
            bigness = 3
        self.power = [0,0]
        medianpower = 70
        minpower = math.floor(medianpower*0.666666)
        maxpower = math.floor(medianpower*1.5)
        self.power[0] = random.randint(minpower,maxpower)*bigness
        self.power[1] = random.randint(minpower,maxpower)*bigness
        if (self.species == "turtle" or self.species == "tortoise" or
            self.species == "roach" or self.species == "snail"):
            self.power[0] *= 0.75
            self.power[1] *= 2
        if (self.species == "dragon" or self.species == "demon" or
            self.species == "angel" or self.species == "kraken" or
            self.species == "hydra"):
            self.power[0] *= 2
            self.power[1] *= 1.5
        if (self.species == "brontosaur" or self.species == "whale" or
            self.species == "albatross" or self.species == "elephant"):
            self.power[0] *= 1.3
            self.power[1] *= 1.6
        self.power[0] = math.floor(self.power[0])
        self.power[1] = math.floor(self.power[1])
        self.measure = "long"
        if self.species in humanoids+bipedalA+bipedalB:
            self.measure = "tall"
        self.measurement = math.floor((self.power[0]+self.power[1])*0.24*random.uniform(0.9,1.1))
        titles = ["devourer","destroyer","primordial beast","ancient",
                  "wanderer","behemoth","titanic beast","colossal beast","goliath",
                  "gargantuan beast","monumental beast","titan","legendary beast"]
        self.title = random.choice(titles) + " "
        e = Event(self.culture.myMap,a=self.age,kind="birth",sub=self,loc=self.location)
        e.importance = math.floor(((self.power[0]/8)+(self.age/8))*random.uniform(0.66666,1.33333))
        self.birthEvent = e
        self.importance += 15
        self.baseImportance = self.importance
    def setPath(self,n):
        self.path = Path(self.location,n)
    def meander(self):
        if random.random() > 0.5:
            return -1
        nm = random.choice(self.location.neighbors)
        j = 0
        if self.terrain == 0:
            while (nm.watery() == 1 or nm.river != None) and j < 100:
                nm = random.choice(self.location.neighbors)
                j = j+1
        if self.terrain == 1:
            while nm.watery() == 0 and j < 100:
                nm = random.choice(self.location.neighbors)
                j = j+1
        self.travel(nm)
    def step(self):
        if self.path == None:
            return -1
        nextNode = self.path.nextNode(self.location)
        if self.profession == "roadbuilder":
            self.location.linkRoads(nextNode)
        self.travel(nextNode)
        if nextNode == self.path.target and self.profession == "roadbuilder":
            self.path = None
            self.erase()
    def travel(self,n):
        if self.location != None and self in self.location.entities:
            self.location.entities.remove(self)
        self.location = n
        self.location.entities.append(self)
    def agePop(self,scl):
        self.stamina = 1
        self.age += scl
        self.importance = clamp(self.importance*1.002,1,(math.log(max(self.baseImportance*1.5,2)))**3)
        if self.dead == 1:
            return -1
        if random.random() > 0.9:
            self.createWork()
        self.offspring()
        if (self.age > self.culture.oldAge or random.random() > 0.99) and (self.immortal == 0 or self.kind == "group"):
            if random.random() < (0.1*self.number) or self.number == 1 :
                roll = self.number*0.9*(random.random())
                if self.number > 1:
                    self.number = math.ceil(roll)
                else:
                    if roll < 0.2:
                        self.number = 0
        if self.number == 0:
            self.die()
        if self.kind == "army":
            self.power[0] = self.number*self.unitbalance[self.profession][0]
            self.power[1] = self.number*self.unitbalance[self.profession][1]
            self.power[0] *= clamp(self.culture.tech["weaponry"]-0.5,1,3)
            self.power[1] *= clamp(self.culture.tech["defense"]-0.5,1,3)
            m = self.culture.value.mainValues
            if "warriors" in m:
                self.power[0] *= 1.1
            if "builders" in m or "metallurgists" in m:
                self.power[1] *= 1.1
            self.power[0] = math.floor(self.power[0])
            self.power[1] = math.floor(self.power[1])
        if self.path != None:
            if self.profession == "roadbuilder":
                self.speed = 1*math.sqrt(self.culture.tech["transportation"])
                if self.path.hasWater() == 1:
                    self.erase()
                    return -1
            numSteps = math.floor(self.speed)
            if random.random() < self.speed-numSteps:
                numSteps += 1
            if self.path.nextNode(self.location) in self.location.roads:
                numSteps += 1
            for eachStep in range(numSteps):
                self.step()
        elif self.kind == "beast":
            self.meander()
        elif self.profession == "explorer" and self.age > 14:
            self.meander()
    def kill(self,n):
        n = math.floor(n)
        self.number = clamp(self.number-n,0,self.number)
    def erase(self):
        self.importance = 0
        self.die(createEvent=False)
        del self.culture.populations[self.name]
    def die(self,createEvent=True):
        if createEvent == True:
            if self.kind == "group" or self.kind == "army":
                deathKind = "disbanding"
                ii = self.importance
            else:
                deathKind = "death"
                ii = ((self.importance*(self.culture.oldAge/self.age))+self.importance)/2
            e = Event(m=self.culture.myMap,a=-1,kind=deathKind,sub=self,loc=self.location)
            e.importance = ii
            self.deathEvent = e
        if self.location != None:
            if self in self.location.entities:
                self.location.entities.remove(self)
        self.dead = 1
        self.deathAge = self.age
    def hasField(self):
        if self.field == None or self.field == "":
            return False
        return True
    def generateFace(self,mode):
        filename = "./generated/face_" + self.justName() + ".gif"
        if self.face == None:
            if len(self.parents) == 0:
                pp = self.culture.generateCultureFace(1)
                res = 192
                self.face = Face(self.culture,p=self,x=res)
                self.face.generateFace1(pp)
            elif len(self.parents) == 1:
                pp = self.parents[0].generateFace(1)
                res = 192
                self.face = Face(self.culture,p=self,x=res)
                self.face.generateFace1(pp)
            else:
                pp1 = self.parents[0].generateFace(1)
                pp2 = self.parents[1].generateFace(1)
                res = 192
                self.face = Face(self.culture,p=self,x=res)
                self.face.generateFace2(pp1,pp2)
            img = Image.new('HSV',(res,res),(255,0,255))
            drawer = ImageDraw.Draw(img)
            self.face.drawSelf(drawer)
            img = img.convert('RGB')
            img.save(filename,"GIF")
        if mode == 0:
            return filename
        else:
            return self.face
    def offspring(self):
        if self.kind != "person":
            return -1
        if self.dead == 1:
            return -1
        if self.location == None:
            return -1
        if self.location.watery():
            return -1
        if self.location.city == None:
            return -1
        primeAge = 40
        ageGap = abs(self.age-primeAge)
        fertility = (2**-(0.5*((ageGap/7)**2))) # From 0 to 1
        if self.culture.society in ["Hegemony","Empire","Imperium","Monarchy"]:
            if len([x for x in self.parents if x == self.culture.leader]) > 0:
                fertility += 0.1
            if self.culture.leader == self and len(self.kids) == 0:
                fertility += 1
        if ageGap > 28:
            fertility = 0
        if random.random()*28 < fertility:
            parentsList = [self]
            newChild = Population(self.culture,a=0,p=1,kind="person",node=self.location,pars=parentsList)
            newChild.age += 1
            newChild.name = (newChild.name[0],self.name[1])
            self.kids.append(newChild)
    def createWork(self):
        if not (self.hasField()):
            return -1
        if self.profession == "roadbuilder":
            return -1
        if self.culture.leader == self and random.random() > 0.25:
            return -1
        if self.number > 1 and random.random() > 0.35:
            return -1
        if self.age < 15:
            return -1
        kind = "book"
        subj = None
        if self.profession in ["artist","philosopher","politician","tactician"]:
            if random.random() > 0.7:
                kind = "story"
        if self.profession in ["artist"]:
            if random.random() > 0.1:
                kind = "piece"
        if self.profession in ["craftsman","blacksmith","engineer"]:
            if random.random() > 0.1:
                kind = random.choice(["weapon","helmet","bodice"])
        if self.profession in ["magician"]:
            if random.random() > 0.3:
                kind = random.choice(["magic"])
        if kind == "magic":
            newSpell = Magic(self)
            return
        t = ""
        if random.random() > 0.35:
            t = random.choice(["event","event","event","event","pop","pop","item","item"])
        # ONLY have subjects of works be from cultures known to the creator.
        # Todo: use the Opinions when they work.
        knownCultures = []
        knownCultures = knownCultures + [x.other for x in self.culture.cultureOpinions.values() if x.knowledge >= 1]
        knownCultures = knownCultures + [self.culture, self.culture, self.culture, self.culture, self.culture]
        knownCulturesNames = [x.name for x in knownCultures]
        if t == "event":
            e = random.choice(self.culture.myMap.events)
            while ((e.subject.culture != self.culture and random.random() < 0.98-(e.importance/400))
                or (e.subject.culture == self.culture and random.random() < 0.75-(e.importance/100))
                or (e.subject.culture.name not in knownCulturesNames)):
                e = random.choice(self.culture.myMap.events)
            subj = e
            e.importance = e.importance*(1+(random.uniform(0.1,0.4)**2))
        if t == "pop":
            p = self.culture.populations[random.choice(list(self.culture.populations.keys()))]
            cc = random.choice([self.culture,self.culture,self.culture,self.culture,self.culture,self.culture,
                                random.choice(knownCultures)])
            while (random.random() < 0.96-(p.importance/100)):
                p = cc.populations[random.choice(list(cc.populations.keys()))]
            subj = p
            p.importance = p.importance*(1+(random.uniform(0.1,0.4)**2))
        if t == "item" and len(self.culture.items) != 0:
            i = random.choice(self.culture.items)
            cc = random.choice([self.culture,self.culture,self.culture,self.culture,self.culture,self.culture,
                                random.choice(knownCultures)])
            while (random.random() < 0.96-(i.importance/100)) and len(cc.items) != 0:
                i = random.choice(cc.items)
            subj = i
            i.importance = i.importance*(1+(random.uniform(0.1,0.4)**2))
        if random.random() > 0.45:
            fff = self.field
        else:
            fff = random.choice(self.culture.myMap.spheres)
        w = Item(k=kind,c=self.culture,f=fff,s=subj,i=self.importance*random.uniform(0.3,1.8),cr=self)
        if w.importance >= 10:
            self.culture.tech[self.field] *= ((random.uniform(0.001,0.01))+1)
        self.works.append(w)
        self.importance = self.importance*random.uniform(1.05,1.2)
    def justName(self):
        s = self.name[0]
        if self.name[1] != "":
            s += " " + self.name[1]
        return s
    def nameFull(self):
        if self.kind == "army":
            return self.name[0]
        if self.title != "":
            s = self.title
            if self.title[-1] != ' ':
                s += " "
            s += self.name[0]
        else:
            if self.profession == None:
                s = self.kind + " " + self.name[0]
            else:
                s = self.profession + " " + self.name[0]
        if self.name[1] != "":
            s += " " + self.name[1]
        return s
    def popNotes(self):
        s = "The "
        s += self.nameFull()
        s += " is a"
        if self.kind == "army":
            s += "n "
        else:
            s += " "
        s += self.kind + ".\n\n"
        s += self.descrip()
        return s
    def descrip(self):
        s = self.pronouns[self.gender].capitalize() + " " + self.toBe[self.gender]
        if self.kind == "group":
            g = "group of, on average, "
        else:
            g = ""
        if self.age < self.culture.myMap.age and self.kind != "army" and self.dead == 0:
            if self.age > 1000:
                s += " a " + g +str(nearestHundred(self.age)) + "-year-old "
            else:
                s += " a " + g +str(self.age) + "-year-old "
        elif self.age < self.culture.myMap.age and self.kind != "army" and self.dead == 1:
            s += " a " + g +str(self.deathAge) + "-year-old "
        elif self.kind != "army":
            s += " an ageless "
        if self.kind == "deity":
            s += "deity of "
            s += synonym(self.associations[0],seedNum(self.name[0]))
            if self.associations[0] == "mountain":
                s += "s"
            s += " from the " + self.culture.name + " culture"
        if self.kind == "person":
            s += "person from the "
            s += self.culture.name + " culture"
        if self.kind == "group":
            s += "people from the "
            s += self.culture.name + " culture.\n"
            s += self.pronouns[self.gender].capitalize()
            s += " is composed of " + str(self.number) + " people"
        if self.kind == "location":
            s += "location of the "
            s += self.culture.name + " culture"
        if self.kind == "beast":
            s += "wild beast named by the "
            s += self.culture.name + " culture.\n"
            s += self.pronouns[self.gender].capitalize() + " is a "
            s += synonym(self.size,seedNum(self.name[0]))
            s += " " + str(self.measurement) + "-meter-" + self.measure
            s += " " + self.species
            s += " and"
            if self.terrain == 0:
                s += " lives only on land"
            elif self.terrain == 1:
                s += " lives only in water"
            else:
                s += " can fly over both land and water"
            s += ".\n"
            s += self.possessive[self.gender].capitalize()
            s += " temperament is " + self.aggression
        if self.kind == "army":
            s += " a " + self.profession + " army"
            s += " composed of " + str(self.number) + " people.\n"
            s += self.profession.capitalize() + " armies are strong against the following army types:\n"
            m = 0
            for j in self.rps[self.profession]:
                if m == 1:
                    s += ", "
                s += j.capitalize()
                m = 1
        s += ".\n"
        if self.profession != None and self.kind != "army":
            if self.number > 1:
                s += "They are "
            else:
                s += self.pronouns[self.gender].capitalize() + " " + self.toBe[self.gender] + " a"
                vows = ['a','e','i','o','u']
                if self.profession[0] in vows:
                    s += "n "
                else:
                    s += " "
            s += self.profession
            if self.number > 1:
                s += "s"
            s += ".\n"
        if len(self.parents) == 0:
            s += ""
        elif len(self.parents) == 1:
            s += self.possessive[self.gender].capitalize() + " parent is"
            s += " the " + self.parents[0].nameFull()
            s += ".\n"
        else:
            s += self.possessive[self.gender].capitalize() + " parents are"
            s += " the " + self.parents[0].nameFull()
            s += " and"
            s += " the " + self.parents[1].nameFull()
            s += ".\n"
        if len(self.kids) == 0:
            s += ""
        elif len(self.kids) == 1:
            s += self.possessive[self.gender].capitalize() + " child is"
            s += " the " + self.kids[0].nameFull()
            s += ".\n"
        else:
            s += self.possessive[self.gender].capitalize() + " children are:\n"
            s += "The " + self.kids[0].nameFull()
            qq = len(self.kids)-1
            for i in range(qq):
                cc = self.kids[i+1]
                if i+1 != qq:
                    s += ",\n the "
                else:
                    s += ",\n and the "
                s += cc.nameFull()
            s += ".\n"
        if self.kind == "beast" or self.kind == "army":
            s += self.possessive[self.gender].capitalize()
            s += " offensive power is approximately equal to "
            s += str(self.power[0]) + " men.\n"
            s += self.possessive[self.gender].capitalize()
            s += " defensive power is approximately equal to "
            s += str(self.power[1]) + " men.\n"
        s += "\n" + self.pronouns[self.gender].capitalize() + " " + self.toBe[self.gender] + " generally considered "
        if self.importance < 13:
            s += "unimportant"
        elif self.importance < 25:
            s += "important"
        elif self.importance < 45:
            s += "very important"
        elif self.importance < 65:
            s += "extremely important"
        else:
            s += "legendary"
        s += ".\n"
        if self.dead == 1:
            s += self.pronouns[self.gender].capitalize()+ " " + self.toBe[self.gender]
            if self.kind == "group":
                s += " disbanded."
            else:
                s += " dead. "
                s += "(Born " + str(self.birthEvent.age) + " years ago; died " + str(self.deathEvent.age) + " years ago)"
        self.description = s
        return s
    def drawHammer(self,drawer,x,y):
        drawer.line([(x-3,y+4),(x+5,y-4)],fill=(0,0,0),width=2)
        p0 = (x+1,y-5)
        p1 = (x+6,y)
        p2 = (x+4,y+2)
        p3 = (x-1,y-3)
        drawer.polygon([p0,p1,p2,p3],fill=self.culture.bannerColor,outline=(0,0,0))
    def drawChevron(self,drawer,x,y,count):
        r = 1
        y = y-math.ceil(count/2)
        for n in range(count):
            if r == 1:
                c = (0,0,0)
            else:
                c = self.culture.bannerColor
            pts = [(x,y-2),(x-3,y+1),(x+3,y+1)]
            drawer.polygon(pts,fill=c)
            y = y+2
            r = 1-r
    def drawSkull(self,drawer,x,y):
        pts = [(x,y),(x-1,y),(x+1,y),(x-2,y),(x+2,y),(x,y+1),
               (x-1,y+1),(x+1,y+1),(x,y-1),(x-2,y-1),(x+2,y-1),
               (x,y-2),(x-1,y-2),(x+1,y-2),(x-2,y-2),(x+2,y-2),(x,y-2)]
        drawer.point(pts,fill=(0,0,255))
        pts = [(x-1,y-1),(x+1,y-1),(x,y+3),(x+2,y+3),(x-2,y+3),
               (x-2,y-3),(x-1,y-3),(x,y-3),(x+1,y-3),(x+2,y-3),
               (x-3,y-2),(x-3,y-1),(x-3,y),(x-3,y+1),(x-2,y+1),
               (x+3,y-2),(x+3,y-1),(x+3,y),(x+3,y+1),(x+2,y+1),
               (x-2,y+2),(x-1,y+2),(x,y+2),(x+1,y+2),(x+2,y+2)]
        drawer.point(pts,fill=(0,0,0))
    def drawShip(self,drawer,x,y):
        c = self.culture.bannerColor
        out = (0,0,0)
        pts = [(x,y),(x,y-6)]
        drawer.line(pts,fill=out,width=1)
        pts = [(x-5,y-5),(x+5,y+5)]
        drawer.chord(pts,start=0,end=180,fill=c,outline=out)
        pts = [(x-3,y-5),(x+3,y-2)]
        drawer.rectangle(pts,outline=out,fill=c)
    def drawPerson(self,drawer,x,y):
        c = self.culture.bannerColor
        out = (0,0,0)
        pts = [(x-3,y-5),(x+3,y+1)]
        drawer.ellipse(pts,outline=out,fill=c)
        pts = [(x-3,y+2),(x+3,y+5)]
        drawer.rectangle(pts,outline=out,fill=c)
    def drawSelf(self,drawer):
        if self.location == None:
            return -1
        if self.location.city != None:
            return -1
        if self.dead == 1:
            return -1
        x = round(self.location.x)
        y = round(self.location.y)
        if self.profession == "roadbuilder":
            self.drawHammer(drawer,x,y)
        elif self.kind == "army":
            self.drawChevron(drawer,x,y,3)
        elif self.kind == "beast":
            self.drawSkull(drawer,x,y)
        elif self.location.watery():
            self.drawShip(drawer,x,y)
        else:
            self.drawPerson(drawer,x,y)
            
        

class Path:
    def __init__(self,n1,n2,land=True):
        self.current = n1
        self.target = n2
        self.nodes = [n1]
        self.path(land)
        self.watery = 0
    def path(self,land):
        a = self.nodes[-1]
        while a != self.target and a != None:
            dist = 10000000
            nn = None
            for p in a.neighbors:
                if p.dist(self.target) < dist and ((p.watery() == 0 and land == True) or land == False):
                    nn = p
                    dist = p.dist(self.target)
            if nn in self.nodes:
                self.watery = 1
                for p in a.neighbors:
                    if p.dist(self.target) < dist:
                        nn = p
                        dist = p.dist(self.target)
            self.nodes.append(nn)
            a = self.nodes[-1]
    def nextNode(self, currentNode):
        if currentNode == None:
            return self.nodes[0]
        nodeIndex = self.nodes.index(currentNode)
        if nodeIndex == len(self.nodes):
            return -1
        return self.nodes[nodeIndex+1]
    def hasWater(self):
        if self.watery == 1:
            return 1
        for n in self.nodes:
            if n.watery() == 1:
                return 1
        return 0
        

class Flag:
    def __init__(self,c):
        self.culture = c
        self.xDim = 384
        self.yDim = 192
        self.colors = []
        self.newColor()
        self.filename = None
    def newColor(self):
        h = random.randint(0,255)
        s = random.randint(32,255)
        v = random.randint(32,255)
        col = (h,s,v)
        self.colors.append(col)
        return col
    def randPt(self):
        pt = (random.randint(0,self.xDim),random.randint(0,self.yDim))
        return pt
    def center(self):
        return (self.xDim/2,self.yDim/2)
    def addTri(self,drawer):
        col = self.newColor()
        p0 = random.choice(self.corners)
        p1 = random.choice([(abs(self.xDim-p0[0]),p0[1]),(p0[0],abs(self.yDim-p0[1]))])
        p2 = random.choice([self.center(),random.choice(self.corners)])
        drawer.polygon([p0,p1,p2],fill=col,outline=col)
    def addRect(self,drawer):
        col = self.newColor()
        p0 = random.choice([(random.randint(0,self.xDim),random.choice([0,self.yDim,self.yDim/2])),
              (random.choice([0,self.xDim,self.xDim/2]),random.randint(0,self.yDim))])
        p1 = random.choice(self.corners)
        drawer.rectangle([p0,p1],fill=col,outline=col)
    def addCirc(self,drawer):
        col = self.newColor()
        p0 = random.choice([(random.randint(0,self.xDim),random.choice([0,self.yDim,self.yDim/2])),
              (random.choice([0,self.xDim,self.xDim/2]),random.randint(0,self.yDim))])
        rad = random.randint(1,self.yDim)
        drawCircle(drawer,p0[0],p0[1],rad,col)
    def getFilename(self):
        if self.filename == None:
            self.genFlag()
        return self.filename
    def genFlag(self):
        img = Image.new('HSV',(self.xDim,self.yDim),self.colors[0])
        drawer = ImageDraw.Draw(img)
        numElements = random.randint(1,3)
        self.corners = [(0,0),(self.xDim,0),(self.xDim,self.yDim),(0,self.yDim)]
        for i in range(numElements):
            element = random.choice(["tri","rect","circ"])
            if element == "tri":
                self.addTri(drawer)
            if element == "rect":
                self.addRect(drawer)
            if element == "circ":
                self.addCirc(drawer)
        drawer.line(self.corners+[(0,0)],fill=(0,0,0),width=8)
        self.filename = "./generated/flag_" + self.culture.name + ".gif"
        img = img.convert('RGB')
        img.save(self.filename,"GIF")

class Language:
    def __init__(self,c):
        self.culture = c
        self.characters()
        self.lengthPref = random.choice([3,5,9])
        self.properNouns = []
        self.name = self.genName()
    def characters(self):
        c = ['b','c','d','f','g','h','j','k','l','m','n','p','q','r','s','t','v','w','x','y','z']
        v = ['a','e','i','o','u']
        e = ['\'','-']
        self.langConsonants = []
        self.langVowels = []
        self.langExtras = []
        count = 32
        n = 1
        while len(c) > 0:
            cons = random.choice(c)
            for l in range(math.floor(count)):
                self.langConsonants.append(cons)
            c.remove(cons)
            n += 1
            count = 32*(1/n)*random.uniform(0.8,1.25)
        count = 32
        while len(v) > 0:
            vow = random.choice(v)
            for l in range(math.floor(count)):
                self.langVowels.append(vow)
            v.remove(vow)
            n += 1
            count = 32*(1/n)*random.uniform(0.8,1.25)
        self.preferredStart = random.choice(self.langConsonants+self.langVowels)
    def isProperNoun(self,word):
        if (word in self.properNouns):
            return True
        for otherCulture in self.culture.myMap.cultures:
            if word in otherCulture.language.properNouns:
                return True
        return False
    def translatePassage(self,passage):
        tokensList = []
        currentToken = ""
        nonWordTokens = ['.',',',' ','?','!','1','2','3','4','5','6','7','8','9','0','\n','\\','-']
        for character in list(passage):
            if character in nonWordTokens:
                if (currentToken != ""):
                    tokensList.append(currentToken)
                    currentToken = ""
                tokensList.append(character)
            else:
                currentToken = currentToken + character
        translatedPassage = ""
        for token in tokensList:
            if token in nonWordTokens:
                translatedPassage = translatedPassage + token
            else:
                addWord = self.translateWord(token)
                if token[0].isupper():
                    addWord = addWord.capitalize()
                    if (self.isProperNoun(token)):
                        addWord = token
                translatedPassage = translatedPassage + addWord
        return translatedPassage
    def translateWord(self,word):
        sd = 0
        ind = 0
        for character in word:
            sd = sd+(ord(character.lower())*getPrime(ind))
        return self.genName(sd)
    def genName(self,sd=None):
        if sd != None:
            random.seed(sd)
        length = random.randint(3,9)
        length = math.floor((length+self.lengthPref)/2)
        n = ""
        con = 0
        vow = 0
        lastchar = '1'
        for k in range(length):
            ctype = random.choice(["con","vow"])
            if vow >= 2:
                ctype = "con"
            if con >= 2:
                ctype = "vow"
            c = lastchar
            if ctype == "con":
                c = random.choice(self.langConsonants)
                vow = 0
                con += 1
            if ctype == "vow":
                c = random.choice(self.langVowels)
                con = 0
                vow += 1
            while c == lastchar and (random.random() > 0.2):
                if ctype == "con":
                    c = random.choice(self.langConsonants)
                    vow = 0
                    con += 1
                if ctype == "vow":
                    c = random.choice(self.langVowels)
                    con = 0
                    vow += 1
            if k == 0 and random.random() < 0.3:
                c = self.preferredStart
            n += c
            lastchar = c
        if sd == None:
            n = n.capitalize()
            self.properNouns.append(n)
        return n

class Map:
    def __init__(self,aAtlas,numNodes,mapDimX,mapDimY):
        self.atlas = aAtlas
        self.n = numNodes
        self.xDim = mapDimX
        self.yDim = mapDimY
        self.landmasses = []
        self.waterBodies = []
        self.regions = []
        self.cities = []
        self.cultures = []
        self.events = []
        self.resourceRegions = []
        self.resourceScale = 1
        self.sealevel = 0.4
        self.date = random.randint(113,974)
        self.seasonStrength = 0.041
        self.setNorth()
        self.biomeColors()
        self.displayNo = None
        self.infoGui = None
        self.extraGui = None
        self.viewmode = 0
        self.drawpops = 1
        self.timeScale = 1
        self.autoCycle = 0
        self.age = random.randint(1000,100000)
        self.year = random.randint(113,974)
        self.roadCol = (20,120,104)
    def relaxAvg(self,strength):
        for i in range(strength):
            for p in self.atlas:
                count = len(p.neighbors)
                xavg = 0
                yavg = 0
                for nn in p.neighbors:
                    xavg += nn.x
                    yavg += nn.y
                xx = xavg/count
                yy = yavg/count
                p.x = xx
                p.y = yy
    def setNorth(self):
        nx = math.floor(random.random()*self.xDim)
        ny = math.floor(random.random()*self.yDim)
        side = random.randint(0,3)
        if side == 0:
            nx = 0
        elif side == 1:
            nx = self.xDim
        elif side == 2:
            ny = 0
        else:
            ny = self.yDim
        self.north = Node(nx,ny,self)
    def nodeLat(self,n):
        return "Latitude: " + str(n.dist(self.north))
    def nodeElevation(self,n):
        return "Elevation: " + str(round((n.elevation-self.sealevel)*self.eScale,1)) + "m"
    def nodeTemp(self,n):
        if n.biome != "water":
            return "Temperature: " + str(round((n.temp*self.tempScale)-30,1)) + " degrees"
        else:
            return "Temperature: " + str(round((n.temp*self.tempScale*0.3),1)) + " degrees"
    def nodeRain(self,n):
        if n.biome != "water":
            rain =  "Rainfall: " + str(round(((n.rainfall)**2)*self.rainfallScale,1)) + "cm/yr"
        else:
            rain = "Rainfall: " + str(round((n.rainfall)*self.rainfallScale*0.03,1)) + "cm/yr"
        return rain
    def nodeBiome(self,n):
        return n.biome
    def nodeRegion(self,n):
        if n.region.culturalNames == {}:
            return "Unnamed " + n.region.biome + "\n"
        else:
            names = ""
            for f in n.region.culturalNames.keys():
                q = n.region.culturalNames[f] + " "
                q += synonym(n.region.biome,seedNum(q))
                q += " (" + f + " language)"
                names += q + "\n"
            return names
    def nodeLandmass(self,n):
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
    def nodeRiver(self,n):
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
    def nodeFertility(self,n):
        return "Soil fertility: "+str(math.floor(n.fertility*self.fertScale))+"%"
    def nodeMetallicity(self,n):
        return "Ground metallicity: "+str(math.floor(n.metallicity*self.metalScale))+"ppm"
    def nodeVegetation(self,n):
        return "Vegetation: " + str(math.floor(n.vegetation*self.vegScale)) + "p/km"+chr(0x00B2)
    def nodeHerbivores(self,n):
        return "Herbivores: " + str(math.floor(n.herbivores*self.wildlifeScale)) + "/km"+chr(0x00B2)
    def nodeCarnivores(self,n):
        return "Carnivores: " + str(math.floor(n.carnivores*self.wildlifeScale)) + "/km"+chr(0x00B2)
    def nodeCityInfo(self,n):
        cityInfo = n.city.cityInfo() + "\n" + "\n"
        cityInfo += n.city.cultureInfo()
        return cityInfo
    def nodeTerritory(self,n):
        territory = ""
        if n.culture != None:
            territory = n.culture.shortName() + " territory"
        return territory
    def nodeResReg(self,n):
        if n.resourceRegion != None:
            reg = "Inside " + n.resourceRegion.rootCity.name + " region" + "\n"
            reg += "Total food available: " + str(math.floor(n.resourceRegion.resources[0]*self.rscScale)) + " t/year \n"
            reg += "Total industrial resources available: " + str(math.floor(n.resourceRegion.resources[1]*self.rscScale)) + " t/year \n"
        return reg
    def infoScales(self):
        self.distScale = 12
        self.eScale = random.randint(2000,3000)
        self.tempScale = 73
        self.rainfallScale = 1756
        self.fertScale = 100
        self.metalScale = 140000
        self.vegScale = 18295
        self.wildlifeScale = 1000
        self.rscScale = 50
    def nodeInfo(self,n):
        self.divWidth = 64
        info = ""
        pops = len(n.entities)
        if pops != 0:
            info += "Number of notable entities at this location: " + str(pops) + "\n"
        if n.city != None:
            info += strDivider(self.divWidth)+"\n"
            info += self.nodeCityInfo(n) + "\n"
        if n.resourceRegion != None:
            info += strDivider(self.divWidth)+"\n"
            info += self.nodeResReg(n)+"\n"
        info += self.nodeTerritory(n) + "\n"
        info += strDivider(self.divWidth)+"\n"
        info += self.nodeRegion(n) + "\n"
        if n.landmass != None:
            info += self.nodeLandmass(n) + "\n"
        if n.river != None:
            info += self.nodeRiver(n) + "\n"
        info += self.nodeElevation(n) + "\n"
        info += self.nodeTemp(n) + "\n"
        info += self.nodeRain(n) + "\n"
        if n.biome != "water":
            info += strDivider(self.divWidth)+"\n"
            info += self.nodeFertility(n) + "\n"
            info += self.nodeMetallicity(n) + "\n"
            info += self.nodeVegetation(n) + "\n"
            info += self.nodeHerbivores(n) + "\n"
            info += self.nodeCarnivores(n) + "\n"
        return info
    def nearestNode(self,xx,yy):
        n = self.atlas[0]
        minDist = 1000000
        search = Node(xx,yy)
        for p in self.atlas:
            dist = p.dist(search)
            if dist < minDist:
                minDist = dist
                n = p
        return n
    def nearestCity(self,xx,yy):
        if len(self.cities) == 0:
            return self.atlas[0]
        n = self.cities[0].node
        minDist = 1000000
        search = Node(xx,yy,self)
        for p in self.cities:
            dist = search.dist(p.node)
            if dist < minDist:
                minDist = dist
                n = p
        return n
    def setSeaLevel(self,n):
        self.sealevel = n
    def perlinElevation(self,octaves,scale=1):
        noiseThing = noiseMaker(math.floor(self.xDim/scale),math.floor(self.yDim/scale))
        for p in self.atlas:
            p.elevation = noiseThing.turbulence(math.floor(p.x),math.floor(p.y),2**octaves)
    def flatten(self):
        for p in self.atlas:
            p.elevation = 0.5
    def elevationAdd(self,amount):
        for p in self.atlas:
            p.elevation = clamp(p.elevation+amount,0,1)
    def randomizeElevation(self,anchor,rng):
        for p in self.atlas:
            p.elevation = anchor-(random.random()*(rng/2))
    def clampElevation(self):
        for p in self.atlas:
            p.elevation = clamp(p.elevation,0,1)
    def smooth(self,strength):
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
                p = clamp(p.elevation,self.sealevel-0.01,255)
    def buildLandmass(self,root):
        if root.landmass != None or root.watery() == 1:
            return -1
        root.landmass = Landmass(root,self)
        self.landmasses.append(root.landmass)
    def buildAllLand(self):
        print("Building landmasses...")
        for p in self.atlas:
            self.buildLandmass(p)
    def buildBodyWater(self,root,waterlev=0):
        if waterlev == 0:
            waterlev = self.sealevel
        if root.landmass != None or root.bodyWater != None:
            return -1
        if root.elevation < waterlev:
            root.bodyWater = bodyWater(root,self.sealevel)
            self.waterBodies.append(root.bodyWater)
    def buildAllWater(self):
        print("Building lakes/oceans...")
        for p in self.atlas:
            self.buildBodyWater(p)
    def buildRegions(self):
        print("Building world regions...")
        for p in self.atlas:
            if p.region == None:
                newReg = Region(p)
                self.regions.append(newReg)
    def addRiver(self,length):
        c = 0
        landmass = self.landmasses[0]
        while c < len(self.landmasses):
            landmass = self.landmasses[math.floor(random.random()*len(self.landmasses))]
            if landmass.size > length*8:
                c = 100000000
            c += 1
        landmass.addRiver(length)
    def addMinorRiver(self,count):
        for i in range(count):
            self.addRiver(self.n/512)
    def addMajorRiver(self,count):
        for i in range(count):
            self.addRiver(self.n/128)
    def cullStreams(self):
        for l in self.landmasses:
            l.cullStreams(self.n/32)
    def addSineHill(self,xx,yy,maximum=0.25,radius=128):
        hillCenter = Node(xx,yy,world)
        for p in self.atlas:
            dist = p.dist(hillCenter)
            if dist <= radius:
                multiplier = distMod(dist,radius)
                if p.elevation > 0.75:
                    multiplier = multiplier*((1-p.elevation)**2)
                p.elevation = p.elevation+(maximum*multiplier)
    def addHill(self,xx,yy,maximum=0.25,radius=128):
        hillCenter = Node(xx,yy,world)
        for p in self.atlas:
            dist = p.dist(hillCenter)
            if dist <= radius:
                multiplier = random.uniform(0.99,1.01)
                p.elevation = maximum*multiplier
    def addMountains(self,num=5,height=0.2):
        avgRad = self.xDim/3.5
        for i in range(num):
            xx = math.floor(random.random()*self.xDim)
            yy = math.floor(random.random()*self.yDim)
            hillRad = avgRad*random.uniform(0.8,1.25)
            hillHeight = height*random.uniform(0.8,1.25)
            self.addSineHill(xx,yy,hillHeight,hillRad)
    def addHills(self,num=8,height=0.1):
        avgRad = self.xDim/6
        for i in range(num):
            xx = math.floor(random.random()*self.xDim)
            yy = math.floor(random.random()*self.yDim)
            hillRad = avgRad*random.uniform(0.8,1.25)
            hillHeight = height*random.uniform(0.8,1.25)
            self.addSineHill(xx,yy,hillHeight,hillRad)
    def addShape(self,shape):
        if shape == "island" or shape == "volcanic":
            self.addSineHill(self.xDim/2,self.yDim/2,0.4,radius=random.uniform(0.6,1.1)*self.xDim/1.5)
            self.smooth(2)
            if shape == "volcanic":
                self.elevationAdd(-0.1)
                self.addSineHill(self.xDim/2,self.yDim/2,0.25,radius=random.uniform(0.6,1.1)*self.xDim*1.5)
                self.smooth(1)
                self.addHill(self.xDim/2,self.yDim/2,random.uniform(0.40,0.49),radius=random.uniform(0.6,1.1)*self.xDim/11)
                self.addSineHill(self.xDim/2,self.yDim/2,-0.2,radius=random.uniform(0.6,1.1)*self.xDim/11)
        if shape == "shore" or shape == "highlands":
            corner = random.randint(0,3)
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
            self.addSineHill(xx,yy,0.45,radius=self.xDim*1.5)
            if shape == "highlands":
                self.addSineHill(xx,yy,0.3,radius=self.xDim)
                self.addMountains()
            self.smooth(2)
        if shape == "archipelago":
            self.addHills(16,0.25)
            self.smooth(2)
        if shape == "plain":
            self.addSineHill(self.xDim/2,self.yDim/2,0.4,radius=self.xDim*5)
            self.smooth(2)
        if shape != "volcanic":
            self.addMountains()
            self.addHills()
        self.smooth(3)
    def addRandomShape(self):
        shp = random.choice(["highlands","plain","volcanic","shore","archipelago","island"])
        self.addShape(shp)
    def erode(self,strength=3):
        for j in range(strength):
            for p in self.atlas:
                if p.elevation < self.sealevel*1.02 and p.elevation > self.sealevel*0.9:
                    if p.hasWaterNeighbor(self.sealevel):
                        p.elevation = p.elevation*0.94
                        for n in p.neighbors:
                            n.elevation = (n.elevation+p.elevation)/2
    def nodeSlopes(self):
        for p in self.atlas:
            p.getSlope()
    def waterdistances(self):
        self.wdmax = 18
        for i in range(self.wdmax):
            for p in self.atlas:
                p.waterdist(self.sealevel)
        for p in self.atlas:
            if p.waterdistance > self.wdmax:
                p.waterdistance = self.wdmax
    def rainfall(self):
        for p in self.atlas:
            rainfall = 0.37*random.uniform(0.9,1.1)
            rainfall = ((rainfall*(0.55-((p.temp*0.72)+(p.elevation*0.3)+(0.8*(p.waterdistance)/(self.wdmax)))))+rainfall)/2
            for n in p.neighbors:
                if n.river != None:
                    rainfall = rainfall*1.2
            if p.river != None:
                rainfall = rainfall*1.2
            for l in p.neighbors:
                if l.watery() == 1:
                    rainfall *= 1.2
            p.rainfall = clamp(rainfall,0,1)
        for p in self.atlas:
            lst = [n.rainfall for n in p.neighbors]
            lst.append(p.rainfall)
            p.rainfall = sum(lst)/len(lst)
    def temperature(self):
        for p in self.atlas:
            s = self.seasonStrength
            seasonMod = clamp(1+(math.sin((math.pi*self.date)/2)*s),1-s,1+s)
            latTemp = 0.6*(p.dist(self.north)/self.xDim)
            elevationTemp = 0.3*(1-p.elevation)
            p.temp = clamp(seasonMod*(elevationTemp+latTemp),0,1)
    def soilProperties(self):
        for p in self.atlas:
            p.metallicity = 0.5*p.elevation*random.uniform(0.666,1.5)
            if p.river != None:
                p.metallicity *= 2
            p.metallicity = clamp(p.metallicity,0,1)
            fertilityBase = (3/clamp(p.waterdistance,1,200)*random.uniform(0.8,1.25))
            p.fertility = fertilityBase
            if p.river != None:
                p.fertility *= 2
            p.fertility = clamp(p.fertility,0,1-(0.25*p.metallicity))
        for p in self.atlas:
            p.metallicity = sum([n.metallicity for n in p.neighbors])/len(p.neighbors)
            p.fertility = sum([n.fertility for n in p.neighbors])/len(p.neighbors)
    def vegetation(self):
        self.temperature()
        self.rainfall()
        for p in self.atlas:
            p.setVegetation()
    def setBiomes(self):
        self.nodeSlopes()
        self.vegetation()
        for p in self.atlas:
            p.setBiome(self.sealevel)
            p.biomeColor = self.biomeColors[p.biome]
            slope = clamp((p.realSlope()*(7500)),-32,32)
            shade = math.floor((-16)+p.biomeColor[2]+slope+(((p.elevation+1)**3)*16))
            p.biomeColor = (p.biomeColor[0],p.biomeColor[1],shade)
    def setWildlife(self):
        for p in self.atlas:
            p.herbivores = clamp((p.vegetation*random.uniform(0.8,1.25))**1.5,0,1)
            p.carnivores = clamp(((p.herbivores*1.2)**1.5),0,1)
    def biomeColors(self):
        bColors = {}
        bColors["desert"] = (16,64,142)
        bColors["savanna"] = (48,64,136)
        bColors["tundra"] = (92,32,104)
        bColors["shrubland"] = (72,96,128)
        bColors["boreal forest"] = (100,64,92)
        bColors["forest"] = (90,108,108)
        bColors["tropical forest"] = (78,162,96)
        bColors["frost"] = (134,32,188)
        bColors["mountain"] = (64,0,96)
        bColors["water"] = Tools.waterColor
        self.biomeColors = bColors
        self.waterHue = 142
        self.waterSaturation = 128
        self.waterValueMin = 44
        self.waterValueMax = 84
    def values(self):
        self.valuesOutputs = {"travelers":0,
                              "craftsmen":0,
                              "traders":0,
                              "superstition":0,
                              "metallurgists":0,
                              "worshippers":0,
                              "freedom":0,
                              "shamans":0,
                              "astrology":0,
                              "materialists":0,
                              "agriculture":0,
                              "collectivists":0,
                              "builders":0,
                              "naturalists":0,
                              "simplicity":0,
                              "greed":0,
                              "sailors":0,
                              "warriors":0}
        self.values = {}
        self.values["swimming"] = {"simplicity":0.2,
                   "sailors":0.5,
                   "astrology":0.15,
                   "freedom":0.4,
                   "warriors":0.1}
        self.values["food"] = {"agriculture":-0.1,
                   "greed":-0.25,
                   "materialists":0.55,
                   "collectivists":0.25,
                   "simplicity":0.45,
                   "worshippers":-0.1,
                   "superstition":-0.2,
                   "traders":0.2,
                   "warriors":0.3,
                   "builders":0.25}
        self.values["darkness"] = {"travelers":0.2,
                   "collectivists":-0.4,
                   "superstition":0.7,
                   "greed":0.45,
                   "astrology":0.25,
                   "materialists":-0.15,
                   "shamans":0.15,
                   "freedom":-0.1,
                   "warriors":0.75,
                   "worshippers":0.5,
                   "builders":0.25}
        self.values["movement"] = {"travelers":0.65,
                   "sailors":0.3,
                   "traders":0.55,
                   "astrology":0.15,
                   "builders":-0.3,
                   "simplicity":0.4,
                   "materialists":-0.15,
                   "freedom":0.85,
                   "naturalists":0.2,
                   "collectivists":-0.2,
                   "warriors":0.3,
                   "greed":0.2}
        self.values["plantlife"] = {"agriculture":0.25,
                   "greed":0.15,
                   "naturalists":0.25,
                   "shamans":0.05,
                   "craftsmen":0.4,
                   "traders":0.15,
                   "builders":0.4,
                   "simplicity":-0.25,
                   "collectivists":0.15}
        self.values["nature"] = {"naturalists":0.45,
                   "shamans":0.2,
                   "agriculture":0.25,
                   "freedom":0.15,
                   "travelers":0.15,
                   "simplicity":0.35,
                   "collectivists":0.2,
                   "superstition":0.3,
                   "astrology":0.1,
                   "metallurgists":0.4,
                   "warriors":0.2,
                   "worshippers":0.2}
        self.values["growth"] = {"shamans":0.1,
                   "agriculture":0.55,
                   "naturalists":0.45,
                   "metallurgists":-0.2,
                   "freedom":0.15,
                   "astrology":-0.2,
                   "collectivists":0.4,
                   "materialists":-0.1,
                   "warriors":0.25,
                   "builders":0.25,
                   "greed":0.1,
                   "worshippers":0.1}
        self.values["sky"] = {"travelers":0.55,
                   "craftsmen":0.3,
                   "traders":0.25,
                   "superstition":0.5,
                   "metallurgists":0.1,
                   "worshippers":0.6,
                   "freedom":0.55,
                   "shamans":0.1,
                   "astrology":0.7,
                   "simplicity":0.4,
                   "builders":0.1}
        self.values["constellations"] = {"astrology":0.6,
                   "superstition":0.3,
                   "worshippers":0.3,
                   "freedom":0.15,
                   "materialists":-0.15,
                   "shamans":0.05}
        self.values["earth"] = {"metallurgists":1.1,
                   "craftsmen":0.9,
                   "traders":0.45,
                   "materialists":0.75,
                   "agriculture":0.25,
                   "collectivists":0.4,
                   "builders":0.65,
                   "freedom":-0.15,
                   "worshippers":-0.1,
                   "greed":0.4,
                   "superstition":-0.2,
                   "warriors":0.4,
                   "astrology":-0.2}
        self.values["fields"] = {"agriculture":0.75,
                   "builders":0.3,
                   "materialists":0.5,
                   "naturalists":0.35,
                   "superstition":-0.2,
                   "simplicity":-0.3,
                   "collectivists":0.25,
                   "warriors":0.3,
                   "freedom":0.1,
                   "astrology":0.1}
        self.values["sunlight"] = {"worshippers":0.9,
                   "astrology":0.6,
                   "naturalists":0.15,
                   "travelers":0.25,
                   "simplicity":0.75,
                   "freedom":0.6,
                   "materialists":-0.2,
                   "warriors":0.5,
                   "builders":-0.1}
        self.values["ice"] = {"superstition":0.3,
                   "simplicity":-0.4,
                   "freedom":-0.4,
                   "travelers":-0.45,
                   "materialists":0.4,
                   "sailors":0.15,
                   "shamans":0.35,
                   "greed":0.35,
                   "metallurgists":0.3,
                   "warriors":0.45,
                   "agriculture":-0.2,
                   "collectivists":0.05}
        self.values["fear"] = {"superstition":0.7,
                   "worshippers":0.5,
                   "shamans":0.55,
                   "freedom":-0.3,
                   "collectivists":0.35,
                   "simplicity":-0.3,
                   "builders":0.4,
                   "greed":0.65,
                   "materialists":-0.25,
                   "warriors":0.9,
                   "astrology":0.2}
        self.values["death"] = {"freedom":0.15,
                   "collectivists":-0.15,
                   "warriors":0.65,
                   "greed":0.45,
                   "travelers":0.1,
                   "builders":-0.1,
                   "simplicity":-0.1,
                   "materialists":0.2}
        self.values["water"] = {"sailors":1.5,
                   "simplicity":0.3,
                   "freedom":0.65,
                   "travelers":0.6,
                   "builders":0.25,
                   "craftsmen":0.1,
                   "astrology":0.4,
                   "traders":0.7,
                   "warriors":0.15}
    def influences(self):
        self.influenceOutputs = {"sky":0,
                                 "sunlight":0,
                                 "fields":0,
                                 "earth":0,
                                 "ice":0,
                                 "fear":0,
                                 "water":0,
                                 "growth":0,
                                 "nature":0,
                                 "plantlife":0,
                                 "movement":0,
                                 "darkness":0,
                                 "swimming":0,
                                 "food":0}
        self.influences = {}
        self.influences["latitude"] = {"sky":0.2,
                       "sunlight":-0.3,
                       "fear":0.2,
                       "movement":0.2,
                       "constellations":0.5}
        self.influences["elevation"] = {"sky":0.3,
                       "sunlight":0.15,
                       "fields":0.1,
                       "earth":0.1,
                       "ice":0.1,
                       "fear":0.05,
                       "movement":0.4,
                       "death":0.1,
                       "constellations":0.1}
        self.influences["hills"] = {"sky":0.1,
                       "sunlight":0.1,
                       "earth":0.1,
                       "fear":0.05,
                       "movement":0.2,
                       "death":0.1,
                       "constellations":0.1}
        self.influences["temperature"] = {"sunlight":0.2,
                       "darkness":-0.1,
                       "sky":0.1,
                       "movement":0.25,
                       "fear":0.15,
                       "water":0.05,
                       "food":0.1,
                       "death":-0.1}
        self.influences["vegetation"] = {"nature":0.15,
                       "food":0.15,
                       "darkness":0.05,
                       "plantlife":0.2,
                       "sunlight":-0.1,
                       "fields":-0.1,
                       "death":0.1}
        self.influences["rainfall"] = {"water":0.4,
                       "growth":0.2,
                       "nature":0.1,
                       "sky":0.25,
                       "movement":0.1,
                       "darkness":0.1,
                       "plantlife":0.1,
                       "swimming":0.15,
                       "death":-0.1}
        self.influences["vegetation"] = {"plantlife":0.6,
                       "growth":0.15,
                       "nature":0.35,
                       "earth":0.2,
                       "fields":0.15}
        self.influences["fertility"] = {"plantlife":0.35,
                       "growth":0.25,
                       "nature":0.15,
                       "earth":0.3,
                       "fields":0.25}
        self.influences["metallicity"] = {"earth":0.7,
                       "darkness":0.15,
                       "nature":0.1,
                       "death":0.1}
        self.influences["carnivores"] = {"darkness":0.35,
                       "food":0.25,
                       "fear":0.4,
                       "nature":0.15,
                       "death":0.45}
        self.influences["herbivores"] = {"nature":0.35,
                       "food":0.65,
                       "growth":0.3,
                       "earth":0.1,
                       "plantlife":-0.1}
        self.influences["rivers"] = {"water":3,
                       "swimming":1.5,
                       "nature":0.1,
                       "food":0.35,
                       "growth":0.1,
                       "movement":0.3}
        self.influences["water"] = {"water":4.5,
                       "swimming":1.5,
                       "sky":0.15,
                       "fear":0.2,
                       "darkness":0.35,
                       "food":0.25,
                       "movement":0.15,
                       "ice":0.1}
        self.influences["forest"] = {"plantlife":0.6,
                       "darkness":0.25,
                       "nature":0.4,
                       "growth":0.2,
                       "earth":0.1,
                       "food":0.3,
                       "fields":0.1,
                       "death":0.05}
        self.influences["desert"] = {"earth":0.4,
                       "sky":0.3,
                       "sunlight":0.7,
                       "food":-0.3,
                       "fields":0.35,
                       "death":0.3}
        self.influences["shrubland"] = {"food":0.25,
                       "fields":0.4,
                       "sky":0.35,
                       "sunlight":0.2,
                       "earth":0.25,
                       "nature":0.1}
        self.influences["savanna"] = {"food":0.2,
                       "fields":0.6,
                       "sky":0.4,
                       "sunlight":0.3,
                       "nature":0.2,
                       "earth":0.1}
        self.influences["tundra"] = {"food":-0.35,
                       "fields":0.35,
                       "fear":0.2,
                       "darkness":0.4,
                       "sky":0.2,
                       "growth":-0.1,
                       "earth":0.3,
                       "ice":0.5,
                       "death":0.2}
        self.influences["mountain"] = {"food":-0.3,
                       "fields":0.1,
                       "sky":0.7,
                       "sunlight":0.35,
                       "earth":0.45,
                       "fear":0.2,
                       "nature":0.15,
                       "ice":0.5,
                       "death":0.1}
        self.influences["tropical forest"] = {"food":0.35,
                       "fear":0.4,
                       "darkness":0.35,
                       "fields":-0.25,
                       "nature":0.5,
                       "plantlife":0.6,
                       "growth":0.45,
                       "earth":0.1,
                       "death":0.25}
        self.influences["boreal forest"] = {"food":0.1,
                       "nature":0.25,
                       "growth":0.1,
                       "fields":0.2,
                       "darkness":0.2,
                       "plantlife":0.25,
                       "earth":0.15,
                       "ice":0.15}
        self.influences["frost"] = {"ice":0.75,
                       "earth":0.2,
                       "food":-0.3,
                       "plantlife":-0.2,
                       "fear":0.25,
                       "darkness":0.2,
                       "fields":0.15,
                       "water":0.1,
                       "death":0.4}
    def technologySetup(self):
        self.technologies = {}
        self.technologies["weaponry"] = 1
        self.technologies["defense"] = 1
        self.technologies["agriculture"] = 1
        self.technologies["production"] = 1
        self.technologies["metallurgy"] = 1
        self.technologies["medicine"] = 1
        self.technologies["government"] = 1
        self.technologies["equality"] = 1
        self.technologies["art"] = 1
        self.technologies["philosophy"] = 1
        self.technologies["research"] = 1
        self.technologies["transportation"] = 1
        self.technologies["magic"] = 1
        self.techtiers = ["primitive","bronze age",
                           "classical period","medieval","pre-industrial",
                           "industrial","contemporary","contemporary","contemporary"]
        # e.g. "medieval humanity", "classical period humanity"
    def biomeWildlife(self):
        self.wildlife = {}
        self.wildlife["desert"] = [["hare","camel"],["scorpion","snake","vulture","lizard"]]
        self.wildlife["savanna"] = [["antelope","elephant","rhinoceros","giraffe"],["lion","hyena"]]
        self.wildlife["tundra"] = [["hare","reindeer","yak"],["wolf","bear"]]
        self.wildlife["shrubland"] = [["hare","deer","tortoise","goat"],["fox","lynx"]]
        self.wildlife["boreal forest"] = [["hare","reindeer","moose","yak"],["wolf","bear","owl","fox"]]
        self.wildlife["forest"] = [["deer","tortoise","tapir","squirrel"],["tiger","bear","lynx"]]
        self.wildlife["tropical forest"] = [["monkey","gorilla","tapir","sloth","frog","hippopotamus"]
        ,["tiger","snake","panther","lizard"]]
        self.wildlife["frost"] = [["penguin","hare"],["bear","wolf"]]
        self.wildlife["mountain"] = [["goat","alpaca","yak"],["wolf","cougar","hawk","snake"]]
        self.wildlife["water"] = [["carp","salmon","tuna","manatee","whale","turtle","lobster","crab"]
        ,["shark","seal","squid","octopus","swordfish","dolphin"]]
    def godSpheres(self):
        s0 = list(self.influences.keys())
        s1 = list(self.values.keys())
        s2 = list(self.valuesOutputs.keys())
        s3 = list(self.technologies.keys())
        self.spheres = s0+s1+s2+s3
    def nearestCityDist(self,xx,yy):
        if len(self.cities) < 1:
            return self.xDim
        else:
            c = self.nearestCity(xx,yy).node
            d = Node(xx,yy,self).dist(c)
            return d
    def placeCity(self,xx,yy,pop=50,culture=None,node=None):
        if node != None:
            cityNode = self.nearestNode(xx,yy)
        else:
            cityNode = node
        if cityNode.biome != "water" and cityNode.city == None:
            newCity = City(cityNode,pop,culture,self)
            return 1
        else:
            return -1
    def randomCity(self):
        cityNode = random.choice(self.atlas)
        tries = 32
        q = 0
        while (cityNode.biome == "water" or cityNode.city != None or
               cityNode.x < 32 or cityNode.x > self.xDim-32 or cityNode.y < 32 
               or cityNode.y > self.yDim-32 or self.nearestCityDist(cityNode.x,cityNode.y) < 32
               or (cityNode.hasWaterNeighbor(self.sealevel) == 0 and q < tries)):
            q = q+1
            cityNode = random.choice(self.atlas)
        newCity = City(cityNode,pop=random.randint(12,136),m=self)
    def scatterCities(self,n):
        for k in self.atlas:
            k.defaultRoads()
        for i in range(n):
            self.randomCity()
    def scatterBeasts(self,n):
        for i in range(n):
            anchor = random.choice(self.atlas)
            while anchor.city != None:
                anchor = random.choice(self.atlas)
            if anchor.watery() == 1:
                aquatic = 1
            else:
                aquatic = 0
            culture = self.nearestCity(anchor.x,anchor.y).culture
            age = random.randint(24,256)
            monster = Population(c=culture,n=None,t="",a=age,p=1,kind="beast",node=anchor)
    def drawGraph(self,gui=None):
        visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawGraph(graphDraw)
        for n in self.atlas:
            n.drawPoint(graphDraw,1,"red")
        visualAtlas.show()
    def drawElevation(self,pts=0,sea=1,gui=None):
        visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawElevation(graphDraw,self.sealevel*sea)
        for n in self.atlas:
            n.drawElevation(graphDraw)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw)
        visualAtlas.show()
    def drawLandmass(self,pts=0,nbrs=0,gui=None):
        visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawLandmass(graphDraw,self.sealevel)
        for n in self.atlas:
            n.drawLandmass(graphDraw,pts,nbrs)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawPath(graphDraw)
        visualAtlas.show()
    def drawWildlife(self,gui=None):
        visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawWildlife(graphDraw,self.sealevel)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw,self.xDim)
        visualAtlas.show()
    def displayNode(self,event):
        clickedNode = self.nearestNode(event.x-2,event.y-2)
        if len(clickedNode.entities) == 0 and clickedNode.city == None:
            for n in clickedNode.neighbors:
                if len(n.entities) > 0:
                    clickedNode = n
                if n.city != None:
                    clickedNode = n
        self.displayString.set(self.nodeInfo(clickedNode))
        self.displayNo = clickedNode
    def redraw(self):
        visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        if self.viewmode == 0:
            rds = []
            for tri in self.triangles:
                tri.drawReal(graphDraw,self.sealevel)
            for n in self.atlas:
                n.drawReal(graphDraw,self.sealevel)
                if len(n.roads) > 0:
                    rds.append(n)
            for l in self.landmasses:
                for r in l.rivers:
                    r.drawRiver(graphDraw,self.xDim)
            for n in rds:
                n.drawRoads(graphDraw,self.roadCol)
            for c in self.cities:
                c.drawSelf(graphDraw)
        elif self.viewmode == 1:
            for tri in self.triangles:
                tri.drawTerritory(graphDraw,self.sealevel)
            for l in self.landmasses:
                for r in l.rivers:
                    r.drawRiver(graphDraw,self.xDim)
            for c in self.cities:
                c.drawSelf(graphDraw)
        if self.drawpops == 1:
            for c in self.cultures:
                c.drawPops(graphDraw)
        visualAtlas = visualAtlas.convert("RGB")
        visualAtlas.save(self.mapname,"GIF")
        photo = Image.open(self.mapname)
        self.img = ImageTk.PhotoImage(photo)
        self.lbl.configure(image = self.img)
        self.lbl.image = self.img
        if self.displayNo != None:
            self.displayString.set(self.nodeInfo(self.displayNo))
    def updateTiming(self):
        self.date += self.timeScale
        self.age += self.timeScale
    def updateResources(self):
        self.vegetation()
        for r in self.resourceRegions:
            r.updateReg()
    def updateTerritory(self):
        for c in self.cities:
            c.node.culture = c.culture
            c.node.allegiance = -1/c.population
        for p in self.atlas:
            p.updateAllegiance(self.sealevel)
        for c in self.cities:
            c.node.culture = c.culture
            c.node.allegiance = -1/c.population
    def updateDemogs(self):
        for l in self.cultures:
            l.updateCulture()
        for c in self.cities:
            c.updateDemog()
    def updatePops(self):
        for c in self.cultures:
            c.updatePops()
    def updateEvents(self):
        for e in self.events:
            e.ageEvent()
    def initOpinions(self):
        for c in self.cultures:
            c.initOpinions()
    def nextTurn(self):
        self.updateTiming()
        self.updateResources()
        self.updateDemogs()
        self.updatePops()
        self.updateEvents()
        self.updateTerritory()
        self.redraw()
    def autoTurnsStart(self):
        self.autoCycle = 1-self.autoCycle
        self.autoTurns()
    def autoTurns(self):
        if self.autoCycle == 1:
            self.nextTurn()
            self.gui.after(1000,self.autoTurns)
    def eventInfo(self,e=None):
        if e == None:
            return -1
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        self.eString = StringVar()
        self.eString.set(e.fullDesc())
        pdsc = Label(self.infoGui,textvariable=self.eString)
        pdsc.pack(anchor=W,side=RIGHT)
        self.displayCulture = e.subject.culture
        self.displayNo = self.displayCulture.origin
        b1 = Button(self.infoGui,text="Society Info",command=self.cultureInfo)
        b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
        c1 = "medium aquamarine"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        g = e.subject
        s = " "+g.justName()+" "
        b1 = Button(self.infoGui,text=s)
        b1.configure(command = lambda self=self, d = g: self.popInfo(d))
        b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
        c1 = "SteelBlue1"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        for g in e.actors:
            s = " "+g.justName()+" "
            b1 = Button(self.infoGui,text=s)
            b1.configure(command = lambda self=self, d = g: self.popInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "SteelBlue2"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def itemInfo(self,e=None):
        if e == None:
            return -1
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        self.eString = StringVar()
        self.eString.set(e.description())
        pdsc = Label(self.infoGui,textvariable=self.eString)
        pdsc.pack(anchor=W,side=RIGHT)
        self.displayCulture = e.culture
        b1 = Button(self.infoGui,text="Society Info",command=self.cultureInfo)
        b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
        c1 = "medium aquamarine"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        g = e.creator
        s = " "+g.justName()+" (Creator) "
        b1 = Button(self.infoGui,text=s)
        b1.configure(command = lambda self=self, d = g: self.popInfo(d))
        b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
        c1 = "SteelBlue1"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        if e.subject != None:
            g = e.subject
            s = " "+g.justName()+" (Subject) "
            b1 = Button(self.infoGui,text=s)
            if g.tt == "pop":
                b1.configure(command = lambda self=self, d = g: self.popInfo(d))
            if g.tt == "event":
                b1.configure(command = lambda self=self, d = g: self.eventInfo(d))
            if g.tt == "item":
                b1.configure(command = lambda self=self, d = g: self.itemInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "SteelBlue2"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def magicInfo(self,e=None):
        if e == None:
            return -1
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        self.eString = StringVar()
        self.eString.set(e.description())
        pdsc = Label(self.infoGui,textvariable=self.eString)
        pdsc.pack(anchor=W,side=RIGHT)
        self.displayCulture = e.culture
        b1 = Button(self.infoGui,text="Society Info",command=self.cultureInfo)
        b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
        c1 = "medium aquamarine"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        g = e.creator
        s = " "+g.justName()+" (Creator) "
        b1 = Button(self.infoGui,text=s)
        b1.configure(command = lambda self=self, d = g: self.popInfo(d))
        b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
        c1 = "SteelBlue1"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def popInfo(self,p=None):
        if p == None:
            return -1
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        self.popString = StringVar()
        self.popString.set(p.popNotes())
        pdsc = Label(self.infoGui,textvariable=self.popString)
        pdsc.pack(anchor=W,side=RIGHT)
        self.displayCulture = p.culture
        if p.kind == "person":
            nn = p.generateFace(0)
            photo = Image.open(nn)
            self.faceImg = ImageTk.PhotoImage(photo)
            self.faceLbl = Label(self.infoGui,image=self.faceImg)
            self.faceLbl.config(borderwidth=32)
            self.faceLbl.photo = photo
            self.faceLbl.pack()
        if p.kind == "city":
            return;
        if p.kind != "beast":
            b1 = Button(self.infoGui,text="Society Info",command=self.cultureInfo)
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "medium aquamarine"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        if p.kind == "deity":
            b1 = Button(self.infoGui,text="Mythology Info",command=self.mythologyInfo)
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "light goldenrod"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        if p.birthEvent != None:
            g = p.birthEvent
            bText = "Birth (" + str(g.age) + " years ago)"
            b1 = Button(self.infoGui,text=bText)
            b1.configure(command = lambda self=self, d = g: self.eventInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "sienna3"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        if p.deathEvent != None:
            g = p.deathEvent
            bText = "Death (" + str(g.age) + " years ago)"
            b1 = Button(self.infoGui,text=bText)
            b1.configure(command = lambda self=self, d = g: self.eventInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "sienna4"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        for g in p.parents:
            s = " "+g.justName()+" (Parent) "
            b1 = Button(self.infoGui,text=s)
            b1.configure(command = lambda self=self, d = g: self.popInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "SteelBlue3"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        for g in p.kids:
            s = " "+g.justName()+" (Child) "
            b1 = Button(self.infoGui,text=s)
            b1.configure(command = lambda self=self, d = g: self.popInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "SteelBlue2"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        for g in p.works:
            s = " "+g.justName()+" (Work) "
            b1 = Button(self.infoGui,text=s)
            b1.configure(command = lambda self=self, d = g: self.itemInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "SteelBlue4"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        for g in p.magic:
            s = " "+g.justName()+" (Magic) "
            b1 = Button(self.infoGui,text=s)
            b1.configure(command = lambda self=self, d = g: self.magicInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "aquamarine"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def cultureInfo(self,reset=False):
        if self.displayNo == None:
            return -1
        if self.displayNo.culture == None:
            return -1
        if reset == True:
            self.displayCulture = self.displayNo.culture
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        photo = Image.open(self.displayCulture.flag.getFilename())
        self.flagImg = ImageTk.PhotoImage(photo)
        self.flagLbl = Label(self.infoGui,image=self.flagImg)
        self.flagLbl.config(borderwidth=16)
        self.flagLbl.photo = photo
        self.flagLbl.pack()
        
        nn = self.displayCulture.generateCultureFace(0)
        photo = Image.open(nn)
        self.faceImg = ImageTk.PhotoImage(photo)
        self.faceLbl = Label(self.infoGui,image=self.faceImg)
        self.faceLbl.config(borderwidth=8)
        self.faceLbl.photo = photo
        self.faceLbl.pack()
        
        self.cultureString = StringVar()
        self.cultureString.set(self.displayCulture.cultureNotes())
        cdsc = Label(self.infoGui,textvariable=self.cultureString)
        cdsc.pack()
        p = self.displayCulture.leader
        s = " The " + p.nameFull() + " "
        b1 = Button(self.infoGui,text="Mythology Info",command=self.mythologyInfo)
        b1.pack(anchor=S,side=BOTTOM,expand=YES,fill=BOTH)
        c1 = "light goldenrod"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        b4 = Button(self.infoGui,text="Technology Info",command=self.techInfo)
        b4.pack(anchor=S,side=BOTTOM,expand=YES,fill=BOTH)
        c1 = "goldenrod3"
        b4.config(bg=c1,activebackground=c1,activeforeground=c1)
        b2 = Button(self.infoGui,text=s)
        b2.configure(command = lambda self=self, e = p: self.popInfo(e))
        b2.pack(anchor=S,side=BOTTOM,expand=YES,fill=BOTH)
        c1 = "SteelBlue2"
        b2.config(bg=c1,activebackground=c1,activeforeground=c1)
        b3 = Button(self.infoGui,text="List all entities of this society",command=self.popListInfo)
        b3.pack(anchor=S,side=BOTTOM,expand=YES,fill=BOTH)
        c1 = "SteelBlue3"
        b3.config(bg=c1,activebackground=c1,activeforeground=c1)
    def popListInfo(self):
        if self.displayCulture == None:
            return -1
        if self.extraGui != None:
            self.extraGui.destroy()
        self.extraGui = Toplevel()
        popsList = [p for p in self.displayCulture.populations.values() if (p.kind not in ["deity","location","beast"])]
        numPops = len(popsList)
        rows = math.ceil(2*math.sqrt(numPops))
        count = 0
        currentFrame = Frame(self.extraGui)
        for i in popsList:
            s = " The "+i.nameFull()+" "
            b1 = Button(currentFrame,text=s)
            b1.configure(command = lambda self=self, d = i: self.popInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=X)
            c1 = "SteelBlue2"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
            count += 1
            if count == rows:
                count = 0
                currentFrame.pack(anchor=S,side=LEFT,expand=YES,fill=Y)
                currentFrame = Frame(self.extraGui)
        currentFrame.pack(anchor=N,side=LEFT,expand=YES,fill=NONE)
    def mythologyInfo(self):
        if self.displayCulture == None:
            return -1
        if self.extraGui != None:
            self.extraGui.destroy()
        self.extraGui = Toplevel()
        self.mythString = StringVar()
        self.mythString.set(self.displayCulture.mythNotes())
        mdsc = Label(self.extraGui,textvariable=self.mythString)
        mdsc.config(justify=LEFT)
        mdsc.pack(anchor=W,side=RIGHT)
        for i in range(len(self.displayCulture.deities)):
            s = " "+self.displayCulture.deities[i].justName()+" "
            b1 = Button(self.extraGui,text=s)
            b1.configure(command = lambda self=self, d = self.displayCulture.deities[i]: self.popInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "SteelBlue2"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def techInfo(self):
        if self.displayCulture == None:
            return -1
        if self.extraGui != None:
            self.extraGui.destroy()
        self.extraGui = Toplevel()
        self.techString = StringVar()
        self.techString.set(self.displayCulture.techNotes())
        mdsc = Label(self.extraGui,textvariable=self.techString)
        mdsc.pack()
    def cityInfo(self):
        if self.displayNo == None:
            return -1
        self.displayCity = self.displayNo
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        self.displayCity.drawTownGen()
        photo = Image.open(self.displayCity.townGen.mapName)
        self.townImg = ImageTk.PhotoImage(photo)
        self.townLbl = Label(self.infoGui,image=self.townImg)
        self.townLbl.config(borderwidth=2)
        self.townLbl.photo = photo
        self.townLbl.pack()
        self.cityString = StringVar()
        structure = self.displayNo.structure()
        if self.displayNo.city != None:
            self.cityString.set(self.displayCity.city.cityNotes())
            cdsc = Label(self.infoGui,textvariable=self.cityString)
            cdsc.pack()
        elif structure != None:
            self.cityString.set(structure)
            cdsc = Label(self.infoGui,textvariable=self.cityString)
            cdsc.pack()
        self.displayCulture = self.displayNo.culture
        b1 = Button(self.infoGui,text="Society Info",command=self.cultureInfo)
        b1.pack(anchor=S,side=BOTTOM,expand=YES,fill=BOTH)
        c1 = "medium aquamarine"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        b2 = Button(self.infoGui,text="Entities Info",command=self.entitiesInfo)
        b2.pack(anchor=S,side=BOTTOM,expand=YES,fill=BOTH)
        c1 = "SteelBlue2"
        b2.config(bg=c1,activebackground=c1,activeforeground=c1)
    def entitiesInfo(self):
        if self.displayNo == None:
            return -1
        entities = self.displayNo.entities
        if entities == []:
            return -1
        if self.extraGui != None:
            self.extraGui.destroy()
        self.extraGui = Toplevel()
        self.popsString = StringVar()
        self.popsString.set("    Notable entities at the selected location:    \n")
        cdsc = Label(self.extraGui,textvariable=self.popsString)
        cdsc.pack()
        for p in entities:
            s = "The " + p.nameFull()
            b1 = Button(self.extraGui,text=s)
            b1.configure(command = lambda self=self, e = p: self.popInfo(e))
            b1.pack(anchor=S,side=BOTTOM,expand=YES,fill=BOTH)
            c1 = "SteelBlue2"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def changeView(self):
        if self.viewmode == 1:
            self.viewmode = 0
        else:
            self.viewmode += 1
        self.redraw()
    def dispEntities(self):
        if self.drawpops == 1:
            self.drawpops = 0
        else:
            self.drawpops += 1
        self.redraw()
    def drawReal(self,gui=None):
        visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        for tri in self.triangles:
            tri.drawReal(graphDraw,self.sealevel)
        for n in self.atlas:
            n.drawReal(graphDraw,self.sealevel)
        for l in self.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw,self.xDim)
        for c in self.cities:
            c.drawSelf(graphDraw)
        if gui == None:
            visualAtlas = visualAtlas.convert("RGB")
            visualAtlas.save("map00.png","PNG")
            visualAtlas.show()
        else:
            self.gui = gui
            self.displayString = StringVar()
            self.displayString.set("No node selected")
            self.infoScales()
            desc = Label(gui,textvariable=self.displayString)
            desc.pack(anchor=W,side=RIGHT)
            visualAtlas = visualAtlas.convert("RGB")
            self.mapname = "./generated/map_" + self.cultures[0].language.genName() + ".gif"
            visualAtlas.save(self.mapname,"GIF")
            photo = Image.open(self.mapname)
            self.img = ImageTk.PhotoImage(photo)
            self.lbl = Label(gui,image=self.img)
            self.lbl.pack(anchor=E,side=LEFT)
            self.lbl.bind("<Button-1>",self.displayNode)
            b0 = Button(gui,text="Next Turn",command=self.nextTurn)
            b0.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "orange red"
            b0.config(bg=c1,activebackground=c1,activeforeground=c1)
            b6 = Button(gui,text=" Toggle auto \n cycle turns ",command=self.autoTurnsStart)
            b6.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "firebrick4"
            b6.config(bg=c1,activebackground=c1,activeforeground=c1)
            b1 = Button(gui,text="Society Info",command=self.cultureInfo)
            b1.configure(command = lambda self=self: self.cultureInfo(reset=True))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "medium aquamarine"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
            b3 = Button(gui,text="Area Info",command=self.cityInfo)
            b3.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "aquamarine"
            b3.config(bg=c1,activebackground=c1,activeforeground=c1)
            b5 = Button(gui,text="Entities Info",command=self.entitiesInfo)
            b5.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "cadet blue"
            b5.config(bg=c1,activebackground=c1,activeforeground=c1)
            b2 = Button(gui,text="Change Mode",command=self.changeView)
            b2.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "salmon2"
            b2.config(bg=c1,activebackground=c1,activeforeground=c1)
            b4 = Button(gui,text="Toggle Entities",command=self.dispEntities)
            b4.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "salmon4"
            b4.config(bg=c1,activebackground=c1,activeforeground=c1)
            self.redraw()
            for i in range(1):
                self.nextTurn()
            gui.mainloop()

#----------------------------------------------------------------------#            
# Let's generate a map

numNodes = 2**14
mapDimX = 960
mapDimY = 960
q = 64
atlas = [Node(-q,-q),Node(mapDimX+q,-q),Node(mapDimY+q,mapDimY+q),Node(-q,mapDimY+q)]
world = Map(atlas,numNodes,mapDimX,mapDimY)

print("Generating points...")
for x in range(numNodes-4):
    nodeX = random.random()*mapDimX
    nodeY = random.random()*mapDimY
    newNode = Node(nodeX,nodeY,world)
    atlas.append(newNode)

npFloatAtlas = np.zeros((numNodes,2))
for q in range(len(atlas)):
    nodeX = atlas[q].x
    nodeY = atlas[q].y
    npFloatAtlas[q] = [nodeX,nodeY]

print("Triangulating...")
from scipy.spatial import Delaunay
triangulation = Delaunay(npFloatAtlas)

trisList = triangulation.vertices
trisVerts = triangulation.points
print("Relaxing points...")
relaxLloyd(npFloatAtlas,2)
for q in range(len(npFloatAtlas)):
    nodeX = npFloatAtlas[q,0]
    nodeY = npFloatAtlas[q,1]
    atlas[q].x = nodeX
    atlas[q].y = nodeY

triangles = []
triIndex = 0
print("Building triangles...")
while triIndex < len(trisList):
    triVertsIndices = trisList[triIndex]
    newTri = Triangle(atlas[triVertsIndices[0]],atlas[triVertsIndices[1]],atlas[triVertsIndices[2]],world)
    triangles.append(newTri)
    triIndex += 1

print("Assigning neighbors...")
for tri in triangles:
    if tri.verts[0].isLinked(tri.verts[1]) == 0:
        tri.verts[0].link(tri.verts[1])
    if tri.verts[1].isLinked(tri.verts[2]) == 0:
        tri.verts[1].link(tri.verts[2])
    if tri.verts[2].isLinked(tri.verts[0]) == 0:
        tri.verts[2].link(tri.verts[0])

world.triangles = triangles

print("Generating terrain...")
world.perlinElevation(6)
world.elevationAdd(-0.35)
world.addRandomShape()
world.setSeaLevel(0.4)
world.cullDots()
world.clampElevation()
world.buildAllWater()
world.buildAllLand()
world.addMajorRiver(10)
world.addMinorRiver(10)
world.waterdistances()
world.soilProperties()
print("Defining biomes...")
world.setBiomes()
world.buildRegions()
world.setWildlife()
world.influences()
world.values()
world.technologySetup()
world.godSpheres()
world.scatterCities(random.randint(8,16))
world.scatterBeasts(random.randint(8,24))
world.initOpinions()
print("Drawing map...")
root = Tk()
world.drawReal(root)