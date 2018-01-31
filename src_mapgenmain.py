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

def synonym(x,seed=0):
    s = {}
    s["mountain"] = ["mountain","peak","ridge"]
    s["savanna"] = ["savanna","plain","prairie"]
    s["shrubland"] = ["shrubland","badlands","bushland"]
    s["forest"] = ["forest","woods","wood","woodland"]
    s["desert"] = ["desert","desert","wastes","barrens"]
    s["tundra"] = ["tundra","steppe","tundra"]
    s["frost tundra"] = ["frost tundra","arctic","alpines","frozen tundra"]
    s["tropical forest"] = ["tropical forest","jungle"]
    s["boreal forest"] = ["boreal forest","woods","wood","taiga"]
    s["carnivores"] = ["carnivores","predators"]
    s["herbivores"] = ["herbivores","livestock"]
    s["fear"] = ["fear","terror"]
    s["warriors"] = ["warriors","fighters","soldiers"]
    s["agriculture"] = ["agriculture","farming","irrigation"]
    s["camp"] = ["bivouac","camp","camp","encampment","campsite"]
    s["village"] = ["village","hamlet"]
    s["township"] = ["township","settlement"]
    s["plantlife"] = ["plantlife","plants","vegetation","flora"]
    s["vegetation"] = ["plantlife","plants","vegetation","flora"]
    s["fields"] = ["fields","farms","pastures"]
    s["metallicity"] = ["metallicity","metals","ore"]
    s["fertility"] = ["fertility","plenty","abundance"]
    s["darkness"] = ["darkness","night","twilight","dusk"]
    s["death"] = ["death","mortality"]
    s["ice"] = ["ice","snow","frost"]
    s["greed"] = ["greed","wealth","gold"]
    s["sky"] = ["sky","stars","heavens"]
    s["collectivists"] = ["collectivists","community","cooperation"]
    s["freedom"] = ["freedom","liberation","liberty"]
    s["book"] = ["book","scroll","volume","document","treatise"]
    s["story"] = ["story","novel","epic","poem"]
    syn = x
    if x in s.keys():
        ch = random.randint(0,len(s[x])-1)
        if seed != 0:
            ch = seed % len(s[x])
        syn = s[x][ch]
    return syn

def seedNum(s):
    v = 0
    for k in s:
        v += ord(k)
    return v

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

def clamp(x,minimum,maximum):
    if x < minimum:
        return minimum
    elif x > maximum:
        return maximum
    else:
        return x

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

class noiseMaker:
    def __init__(self,w,h):
        self.noise = np.random.rand(w,h)
        self.width = w
        self.height = h
    def smoothNoise(self,xx,yy):
        fracX = xx % 1
        fracY = yy % 1
        x1 = (math.floor(xx)+self.width) % self.width
        y1 = (math.floor(yy)+self.height) % self.height
        x2 = (x1+self.width-1) % self.width
        y2 = (y1+self.height-1) % self.height
        tileVal = 0
        tileVal += fracX*fracY*self.noise[x1,y1]
        tileVal += (1-fracX)*fracY*self.noise[x2,y1]
        tileVal += fracX*(1-fracY)*self.noise[x1,y2]
        tileVal += (1-fracX)*(1-fracY)*self.noise[x2,y2]
        return tileVal
    def turbulence(self,xx,yy,size):
        tileVal = 0
        initialSize = size
        while size >= 1:
            tileVal += self.smoothNoise(xx/size,yy/size)*size
            size = size/2
        return tileVal/(initialSize*2)
    
class Node:
    def __init__(self,xx,yy):
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
        self.key = 0
        self.roads = []
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
        self.neighbors.append(newNeighbor)
        newNeighbor.neighbors.append(self)
    def getKey(self):
        return self.key
    def watery(self,sealevel):
        if self.elevation < sealevel or self.river != None:
            return 1
        else:
            return 0
    def hasWaterNeighbor(self,sealevel):
        for n in self.neighbors:
            if n.elevation < sealevel:
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
    def smooth(self):
        nbrs = []
        nbrs.append(self.elevation)
        for i in self.neighbors:
            nbrs.append(i.elevation)
        self.elevation = sum(nbrs)/len(nbrs)
    def defaultRoads(self):
        self.roads = []
        for k in self.neighbors:
            self.roads.append(k)
    def setVegetation(self):
        tempFitness = 1-abs(0.3-self.temp)
        elevationFitness = 1-abs(0.45-self.elevation)
        fertilityFitness = self.fertility+(self.metallicity*0.25)
        rainFitness = 1-abs(self.rainfall-0.8)
        vegFitness = ((tempFitness+1)*(elevationFitness+0.5)*(fertilityFitness+0.5)*(rainFitness+2))
        self.vegetation = clamp(vegFitness/16,0,1)
    def setBiome(self,sl):
        if self.elevation > 0.864:
                self.biome = "mountain"
        elif self.temp < 0.15:
            if self.rainfall < 0.08:
                self.biome = "frost"
            elif self.rainfall < 0.3:
                self.biome = "tundra"
            else:
                self.biome = "shrubland"
        elif self.temp < 0.3:
            if self.rainfall < 0.03:
                self.biome = "tundra"
                if self.temp > 0.5:
                    self.biome = "desert"
            elif self.rainfall < 0.06:
                self.biome = "shrubland"
            elif self.rainfall < 0.3:
                self.biome = "boreal forest"
            elif self.rainfall < 0.4:
                self.biome = "forest"
            else:
                self.biome = "tropical forest"
        elif self.temp < 0.55:
            if self.rainfall < 0.03:
                self.biome = "desert"
            elif self.rainfall < 0.06:
                self.biome = "savanna"
            elif self.rainfall < 0.09:
                self.biome = "shrubland"
            elif self.rainfall < 0.2:
                self.biome = "forest"
            else:
                self.biome = "tropical forest"
        else:
            if self.rainfall < 0.03:
                self.biome = "desert"
            elif self.rainfall < 0.06:
                self.biome = "savanna"
            else:
                self.biome = "tropical forest"
        if self.elevation < sl:
            self.biome = "water"
    def claim(self,n,sealevel=0.4):
        n.culture = self.culture
        inc = 1
        if n.watery(sealevel) == 1:
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
                roll = random.random()
                if roll <= chance:
                    self.claim(n,sealevel)
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
        

class Triangle:
    def __init__(self,aa,bb,cc):
        self.verts = [None,None,None]
        self.verts[0] = aa
        self.verts[1] = bb
        self.verts[2] = cc
        self.neighbors = []
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
            if f.elevation < sealevel:
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
            if f.elevation < sealevel:
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
            if f.elevation < sealevel:
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
            if f.elevation < sealevel:
                underwater = 1
        if underwater == 1:
            v = 32
        dCol = (h,s,v)
        drawer.polygon([self.verts[0].coords(),self.verts[1].coords(),self.verts[2].coords()],fill=dCol,outline=dCol)
    def drawReal(self,drawer,sl):
        elevationList = [self.verts[f].elevation for f in range(len(self.verts))]
        elevation = sum(elevationList)/3
        col = math.floor(elevation*255)
        underwater = 0
        avgHue = math.floor(sum([self.verts[f].biomeColor[0] for f in range(len(self.verts))])/3)
        avgSat = math.floor(sum([self.verts[f].biomeColor[1] for f in range(len(self.verts))])/3)
        avgVal = math.floor(sum([self.verts[f].biomeColor[2] for f in range(len(self.verts))])/3)
        dCol = (avgHue,avgSat,avgVal)
        for f in self.verts:
            if f.elevation < sl:
                underwater = 1
        if underwater == 1:
            dCol = (142,128,clamp(col,64,76))
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
            for m in current.neighbors:
                if (m.elevation > choice.elevation
                    and m.river == None and m not in self.nodes):
                    choice = m
            for k in current.neighbors:
                if (k.elevation < choice.elevation and k.elevation >= current.elevation 
                    and k.river == None and k not in self.nodes):
                    choice = k
            if current == choice:
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
        dCol = (142,128,76)
        for i in range(len(self.nodes)-1):
            n = self.nodes[i]
            n1 = self.nodes[i+1]
            scale = xDim/2
            w = clamp(1/n.slope,0,2048)/scale
            w1 = clamp(1/n1.slope,0,2048)/scale
            drawCircle(drawer,n.x,n.y,w,dCol)
            drawCircle(drawer,n1.x,n1.y,w1,dCol)
            drawTrapezoid(drawer,n.x,n.y,n1.x,n1.y,w,w1,dCol)

class bodyWater:
    def __init__(self,rootNode,sLevel):
        self.sealevel = sLevel
        self.root = rootNode
        self.nodes = []
        self.addNode(rootNode)
        self.fill()
        self.culturalNames = {}
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
    def __init__(self,rootNode,sLevel):
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
    def centermass(self):
        xTotal = sum([p.x for p in self.nodes])
        yTotal = sum([p.y for p in self.nodes])
        xx = xTotal/self.size
        yy = yTotal/self.size
        self.centroid = Node(xx,yy)
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
                if k.elevation < self.sealevel:
                    self.addBoundary(p)
                while k not in self.nodes and k.elevation >= self.sealevel:
                    self.addNode(k)
    def addRiver(self,length):
        root = self.boundary[math.floor(random.random()*(len(self.boundary)-1))]
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
        for o in self.myMap.influenceOutputs.keys():
            self.influenceOutput[o] = 0
        self.translateInfluence()
    def setOutput(self,influenceType,strength):
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
        strength = self.node.herbivores
        inf = "carnivores"
        self.setOutput(inf,strength)
        strength = self.node.carnivores
        inf = "herbivores"
        self.setOutput(inf,strength)
        strength = self.node.temp
        inf = "temperature"
        self.setOutput(inf,strength)
        strength = self.node.elevation
        inf = "elevation"
        self.setOutput(inf,strength)
        strength = self.node.slope
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
        return mVals
    def translateValues(self):
        for q in self.influences.influenceOutput.keys():
            modifier = self.influences.influenceOutput[q]
            roll = random.random()
            if roll > 0.99:
                modifier = 8
            elif roll < 0.01:
                modifier = 0.01
            for v in self.myMap.values[q].keys():
                self.valuesOutput[v] += self.myMap.values[q][v]*modifier*random.uniform(0.8,1.25)

class City:
    def __init__(self,n,pop=50,cltr=None,m=None):
        self.myMap = m
        self.myMap.cities.append(self)
        self.node = n
        self.node.city = self
        self.population = pop
        self.culture = None
        if cltr == None:
            self.culture = Culture(self.node,m)
        else:
            self.culture = cltr
        self.cityType = self.cType(self.population)
        self.name = self.culture.language.genName()
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
        self.age = 0
        self.roads = []
    def updateDemog(self):
        rscShare = self.population/self.node.resourceRegion.totalPop
        rscMax = [0,0]
        rscMax[0] = rscShare*self.node.resourceRegion.resources[0]
        rscMax[1] = rscShare*self.node.resourceRegion.resources[1]
        mpo = random.uniform(0.097,0.101)   # Maximum personal output (max resources production per person)
        m = self.culture.value.mainValues
        if "agriculture" in m:
            mpo *= 1.1
        if "simplicity" in m:
            mpo *= 0.85
        if "warriors" in m:
            mpo *= 1.05
        if "builders" in m:
            mpo *= 1.05
        if "metallurgists" in m:
            mpo *= 1.05
        if "craftsmen" in m:
            mpo *= 1.05
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
        self.foodConsumption = mpc*self.population
        diff = self.foodProduction-self.foodConsumption
        growth = clamp(diff/mpc,-self.population*0.1,self.population*0.1)
        if self.population < 20:
            growth = clamp(growth,1,100)
        self.population = math.ceil(self.population*0.99)  # Age related death
        self.population = clamp(math.floor(self.population+growth+random.choice([-1,0,0,0,0,1])),1,1000000)
        self.cityType = self.cType(self.population)
        self.diaspora()
        for r in self.roads:
            r.build()
    def diaspora(self):
        self.age += 1
        self.threshold = 1000*(0.9996**self.age)
        roll = random.random()
        minRoll = 0.6
        superRoll = 0.99
        m = self.culture.value.mainValues
        if "simplicity" in m:
            self.threshold *= 0.4
            roll = roll*0.5
        if "travelers" in m:
            self.threshold *= 0.8
            roll = roll*1.5
        if "freedom" in m:
            self.threshold *= 0.8
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
            self.threshold *= 1.25
        if (((self.population > self.threshold and roll > minRoll) or roll > superRoll)
        and self.age > 12):
            self.migrate()
    def migrate(self):
        emigrants = math.ceil(self.population*random.uniform(0.1,0.5))
        rng = 64
        xx = clamp(self.node.x + random.randint(-rng,rng),16,self.myMap.xDim)
        yy = clamp(self.node.y + random.randint(-rng,rng),16,self.myMap.yDim)
        n = self.myMap.nearestNode(xx,yy)
        if ((n.culture == self.culture or n.culture == None) 
            and n.city == None and n.resourceRegion == None and n.landmass != None):
            cc = City(n,pop=emigrants,cltr=self.culture,m=self.myMap)
            rd = Road(self.node,n,self)
            self.population = clamp(self.population-emigrants,1,10000000)
            self.age = 0
    def cType(self,p):
        if p <= 35:
            c = synonym("camp")
        elif p <= 200:
            c = synonym("village")
        elif p <= 1000:
            c = synonym("township")
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
        p0 = (xx-5,yy-3)
        p1 = (xx-5,yy+1)
        p2 = (xx-3,yy+1)
        p3 = (xx-3,yy-5)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx+5,yy-4)
        p1 = (xx+5,yy+1)
        p2 = (xx+3,yy+1)
        p3 = (xx+3,yy-4)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx-1,yy-7)
        p1 = (xx-1,yy+2)
        p2 = (xx+1,yy+2)
        p3 = (xx+1,yy-7)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx-6,yy+2)
        p1 = (xx-6,yy-1)
        p2 = (xx+6,yy-1)
        p3 = (xx+6,yy+2)
        drawer.polygon([p0,p1,p2,p3],outline=out,fill=col)
        p0 = (xx,yy)
        p1 = (xx,yy+1)
        drawer.line([p1,p0],fill=out,width=1)
    def drawSelf(self,drawer):
        for r in self.roads:
            r.drawSelf(drawer)
        col = self.culture.bannerColor
        out = (0,0,0)
        if self.population <= 35:
            self.drawTent(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 200:
            self.drawHut(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 1000:
            self.drawVillage(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 5000:
            self.drawTown(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 20000:
            self.drawCity(drawer,self.node.x,self.node.y,col,out)
    def drawTownGen(self):
        self.townGen = Town(self.node,self.myMap,self.name)
        townImg = Image.new("HSV",(self.townGen.xDim,self.townGen.yDim),(0,0,255))
        graphDraw = ImageDraw.Draw(townImg)
        self.townGen.drawSelf(graphDraw)
        townImg = townImg.convert("RGB")
        townImg.save(self.townGen.mapName,"GIF")
    def cityNotes(self):
        s = self.name + "\n\n"
        s += self.cType(self.population) + "\n\n"
        s += "Population: " + str(self.population) + "\n\n"
        return s

class Road:
    def __init__(self,n0,n1,c):
        self.city = c
        self.city.roads.append(self)
        self.start = n0
        self.destination = n1
        self.current = n0
        self.nodes = []
        self.nodes.append(self.start)
        self.size = 2
        self.built = 0
    def build(self):
        if self.built == 1:
            return -1
        rate = 1
        i = 0
        while i < rate and self.current != self.destination:
            i += 1
            n = random.choice(self.current.neighbors)
            while n in self.nodes or n.elevation < self.city.myMap.sealevel:
                n = random.choice(self.current.neighbors)
            for p in self.current.neighbors:
                if p.dist(self.destination) < n.dist(self.destination):
                    if p.elevation > self.city.myMap.sealevel:
                        n = p
                if p == self.current:
                    n = p
            self.linkRoad(self.current,n)
            self.current = n
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
        dCol = (16,128,64)
        nds = [x.coords() for x in self.nodes]
        drawer.line(nds,fill=dCol,width=self.size)

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
                rawAnimal += 0.1
                rawPlant += 0.1
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
        scale = self.myMap.resourceScale
        self.resources[0] = scale*(rawPlant + (rawAnimal*0.7))
        self.resources[1] = scale*(rawMetal + (rawAnimal*0.3) + (rawPlant*0.4))
        if "craftsmen" in m:
            self.resources[1] *= 1.1
        if "builders" in m:
            self.resources[1] *= 1.1
        if "collectivists" in m:
            self.resources[0] *= 1.05
    def updateReg(self):
        for p in self.nodes:
            if p.city != None:
                p.resourceDist = math.sqrt(p.city.population)
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
        self.language = Language(self)
        self.name = self.language.name
        self.populations = {}
        self.pantheon = []
        self.generateMythology()
        self.title = self.setTitle()
        self.flag = Flag(self)
        self.bannerColor = self.flag.colors[0]
        self.leaderTitle = self.setLeaderTitle()
        self.leader = Population(self,t=self.leaderTitle,p=self.leaderCount)
        self.totalPop = self.populationCount()
        self.oldAge = 80
    def populationCount(self):
        p = 0
        for f in self.myMap.cities:
            if f.culture == self:
                p += f.population
        return p
    def updatePops(self):
        for f in self.populations.keys():
            self.populations[f].agePop(self.myMap.timeScale)
        for d in self.deities:
            d.agePop(self.myMap.timeScale)
    def setSociety(self):
        m = self.value.mainValues
        if "warriors" in m and "collectivists" in m and "worshippers" in m:
            return "Hegemony"
        if "greed" in m and "worshippers" in m and "builders" in m:
            return "Empire"
        if "warriors" in m and "greed" in m and "builders" in m and ("travelers" in m or "sailors" in m):
            return "Imperium"
        if "builders" in m and "collectivists" in m and "materialists" in m:
            return "Socialists"
        if "travelers" in m and "sailors" in m and "traders" in m:
            return "Traders"
        if "freedom" in m and "collectivists" in m and "simplicity" in m:
            return "Paleolithic tribe"
        if ("travelers" in m or "sailors" in m) and ("traders" in m or "greed" in m) and "freedom" in m:
            return "Independent merchants"
        if "collectivists" in m and "agriculture" in m and "materialists" in m:
            return "Agricultural communists"
        if "collectivists" in m and "agriculture" in m and "simplicity" in m:
            return "Farming commune"
        if "worshippers" in m and "warriors" in m and ("superstition" in m or "collectivists" in m):
            return "Religious zealots"
        if "shamans" in m and ("naturalists" in m or "astrology" in m) and "superstition" in m:
            return "Shamanic tribe"
        if "shamans" in m and "warriors" in m and ("astrology" in m or "superstition" in m or "worshippers" in m):
            return "Shamanistic warriors"
        if "metallurgists" in m and "builders" in m and "craftsmen" in m and "materialists" in m:
            return "Cooperative artisans"
        if "metallurgists" in m and "builders" in m and "craftsmen" in m and "traders" in m:
            return "Merchant artisans"
        if "freedom" in m and "greed" in m and "traders" in m:
            return "Liberal capitalists"
        if "freedom" in m and "traders" in m and ("builders" in m or "craftsmen" in m or "metallurgists" in m):
            return "Liberal merchant-artisans"
        if "builders" in m and "agriculture" in m and "traders" in m:
            return "Mercantile folk"
        if "builders" in m and "agriculture" in m and ("travelers" in m or "sailors" in m):
            return "Township builders"
        if ("craftsmen" in m or "metallurgy" in m or "builders" in m) and "simplicity" in m and ("naturalists" in m or "shamans" in m):
            return "Naturalist artisans"
        if ("craftsmen" in m or "metallurgy" in m or "builders" in m) and "simplicity" in m and ("superstition" in m or "astrology" in m):
            return "Traditionalist artisans"
        if "greed" in m and "sailors" in m and "warriors" in m:
            return "Pirates"
        if (("travelers" in m or "sailors" in m) and "greed" in m and "warriors" in m):
            return "Raiders"
        if (("travelers" in m or "sailors" in m) and "greed" in m and "simplicity" in m):
            return "Scavengers"
        if "travelers" in m and "simplicity" in m and "freedom" in m:
            return "Hunter-gatherer tribe"
        if "astrology" in m and "superstition" in m and "worshippers" in m:
            return "Religious sovereignty"
        if "collectivists" in m and "agriculture" in m and "naturalists" in m:
            return "Agriculturalists"
        if "travelers" in m and "simplicity" in m and ("metallurgists" in m or "craftsmen" in m or "builders" in m):
            return "Nomadic artisans"
        if "travelers" in m and "simplicity" and ("naturalists" in m or "superstition" in m or "astrology" in m or "shamans" in m):
            return "Nomadic tribe"
        if "materialists" in m and "metallurgists" in m and ("craftsmen" in m or "builders" in m):
            return "Blacksmiths"
        if "warriors" in m and "collectivists" in m and ("travelers" in m or "sailors" in m):
            return "Revolutionary commune"
        if "builders" in m and "metallurgists" in m and "craftsmen" in m and "agriculture" in m:
            return "Syndicalists"
        if "warriors" in m and "builders" in m and ("worshippers" in m or "superstition" in m) and "freedom" not in m and "traders" not in m:
            return "Nationalists"
        if "collectivists" in m and "freedom" in m:
            return "Social Liberals"
        if (("astrology" in m and "superstition" in m) or 
            ("superstition" in m and "worshippers" in m) or
            ("woshippers" in m and "astrology" in m)):
            return "Religious sovereignty"
        if "traders" in m and "freedom" in m:
            return "Merchants"
        if "traders" in m and "collectivists" in m:
            return "Co-operative"
        if "traders" in m and "greed" in m:
            return "Capitalists"
        if "collectivists" in m:
            return "Communalists"
        if ("freedom" in m and "materialists" in m and ("superstition" not in m)):
            return "Scholars"
        if ("materialists" in m and "astrology" in m and 
            ("superstition" not in m or "worshippers" not in m)):
            return "Astronomers"
        if "shamans" in m:
            return "Shamans"
        if "simplicity" in m:
            return "Tribe"
        if "agriculture" in m:
            return "Agriculturalists"
        if "worshippers" in m:
            return "Religious sovereignty"
        if "freedom" in m:
            return "Liberals"
        return "Mixed society"
    def setTitle(self):
        if self.society == "Nationalists":
            return "Nation"
        if self.society == "Religious sovereignty" or self.society == "Religious zealots":
            return random.choice(["Theocracy","Ecclesiarchy","Order"])
        if self.society == "Agriculturalists" or self.society == "Farming commune" or self.society == "Agricultural communists":
            return random.choice(["Farmers","Yeomen","Peasants"])
        if self.society == "Imperium" or self.society == "Hegemony" or self.society == "Empire":
            return "Empire"
        if self.society == "Nomadic artisans" or self.society == "Nomadic peoples" or self.society == "Scavengers":
            return "Nomads"
        if (self.society == "Liberal capitalists" or self.society == "Liberal merchant-artisans" or self.society == "Merchant artisans" or 
            self.society == "Traders" or self.society == "Independent merchants" or self.society == "Mercantile folk"
            or self.society == "Merchants" or self.society == "Capitalists"):
            return random.choice(["Caravans","Proprietors","Holdings"])
        if (self.society == "Blacksmiths" or self.society == "Traditionalist artisans" or
            self.society == "Naturalist artisans" or self.society == "Cooperative artisans"):
            return "Artisans"
        if (self.society == "Socialists" or self.society == "Syndicalists" or self.society == "Revolutionary commune"
            or self.society == "Communalists" or self.society == "Co-operative"):
            return random.choice(["People's Union","Union","Collective"])
        if (self.society == "Shamanistic warriors" or self.society == "Shamanic tribe"
            or self.society == "Shamans"):
            return "Mystics"
        if self.society == "Pirates" or self.society == "Raiders":
            return "Brigands"
        if self.society == "Social Liberals" or self.society == "Liberals":
            return "Republic"
        if self.society == "Scholars" or self.society == "Astronomers":
            return random.choice(["Institute","Academy","College"])
        return "People"
    def setLeaderTitle(self):
        self.leaderCount = 1
        s2 = ""
        if self.society == "Nationalists":
            return (random.choice(["Master","High","Lord",""])+ " " +
                    random.choice(["Commissioner","Chancellor","Harbinger"]))
        if self.society == "Religious sovereignty" or self.society == "Religious zealots":
            return (random.choice(["Grand","High","Supreme","Holy"]) + " " +
                    random.choice(["Pontiff","Priest","Shepherd"]))
        if self.society == "Agriculturalists" or self.society == "Farming commune" or self.society == "Agricultural communists":
            self.leaderCount = random.choice([1,random.randint(2,10)])
            if self.leaderCount == 1:
                s = ""
            else:
                s = random.choice(["Council","Assembly","Soviet","Conference","Directorate"]) + " of "
                s2 = "s"
            return (random.choice(["Head","Chief","Master"])+ " " +
                    random.choice(["Farmer","Agronomist","Foreman"])+s2)
        if self.society == "Imperium" or self.society == "Hegemony" or self.society == "Empire":
            return "Emperor"
        if self.society == "Nomadic artisans" or self.society == "Nomadic peoples" or self.society == "Scavengers":
            return (random.choice(["Chief","Head","Elder"])+ " " +
                    random.choice(["Captain","Dignitary","Herald"]))
        if (self.society == "Liberal capitalists" or self.society == "Liberal merchant-artisans" or self.society == "Merchant artisans" or 
            self.society == "Traders" or self.society == "Independent merchants" or self.society == "Mercantile folk"
            or self.society == "Merchants" or self.society == "Capitalists"):
            self.leaderCount = random.choice([1,random.randint(2,10)])
            if self.leaderCount == 1:
                s = ""
            else:
                s = random.choice(["Cabinet","Assembly","Board","Committee"]) + " of "
                s2 = "s"
            return (s + random.choice(["Primary","Head","Chief",""])+ " " +
                    random.choice(["Executive","Director","Superintendent"]) + s2)
        if (self.society == "Blacksmiths" or self.society == "Traditionalist artisans" or
            self.society == "Naturalist artisans" or self.society == "Cooperative artisans"):
            self.leaderCount = random.choice([1,random.randint(2,10)])
            if self.leaderCount == 1:
                s = ""
            else:
                s = random.choice(["Council","Assembly","Congress"]) + " of "
                s2 = "s"
            return (s + random.choice(["Master","Elder","Grandmaster",""])+ " " +
                    random.choice(["Artificer","Builder","Craftsperson"]) + s2)
        if (self.society == "Socialists" or self.society == "Syndicalists" or self.society == "Revolutionary commune"
            or self.society == "Communalists" or self.society == "Cooperative" or self.society == "Scholars"):
            self.leaderCount = random.choice([1,random.randint(2,20)])
            if self.leaderCount == 1:
                s = ""
            else:
                s = random.choice(["Council","Assembly","Soviet","Conference","Directorate"]) + " of "
                s2 = "s"
            return (s+random.choice(["Prime","Chief","Central",""])+ " " +
                    random.choice(["Director","Governer","Speaker","Chairperson"])+s2)
        if (self.society == "Shamanistic warriors" or self.society == "Shamanic tribe"
            or self.society == "Shamans"):
            return (random.choice(["Elder","High","Grand","Ancestral"])+ " " +
                    random.choice(["Medicine Man","Seer","Shaman"]))
        if self.society == "Pirates" or self.society == "Raiders":
            return (random.choice(["Chief","Head",""])+ " " +
                    random.choice(["Captain","Commander",""]))
        if self.society == "Social Liberals" or self.society == "Liberals":
            self.leaderCount = random.choice([1,1,random.randint(2,10)])
            if self.leaderCount == 1:
                s = ""
            else:
                s = random.choice(["Congress","Chamber","Parliament","Ministry"]) + " of "
                s2 = "s"
            return random.choice(["President","Speaker","Minister"]) + s2
        return "Chief"
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
        self.mythAge = ageOf
        rr = 1.2
        f = Population(self,n=self.language.genName(),t=titles[0],a=ageOf,kind="location")
        f.gender = 0
        if ageOf < self.myMap.age:
            e = Event(self.myMap,ageOf,"genesis",f,[])
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
                nom = self.language.genName()
                if random.random() < (0.25 + (k*0.1)):
                    nom += " " + self.language.genName()
                ent = Population(self,n=nom,t=titles[k],a=age,kind="deity")
                ent.kids = []
                for u in range(parent):
                    pp = random.choice(self.deities)
                    while pp in ent.parents:
                        pp = random.choice(self.deities)
                    ent.parents.append(pp)
                    pp.kids.append(ent)
                ent.associate()
                self.deities.append(ent)
                e = Event(self.myMap,age,"birth",ent,ent.parents)
                e.importance = 100/k
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
        """s = self.name + " Mythology" + "\n"
        s += strDivider(128)"""
        for d in self.deities:
            if self.deities.index(d) != 0:
                s += "\n\n"
            s += "[" + str(self.deities.index(d)) + "]     "
            if d.age > self.myMap.age:
                s += "Since the beginning of time, "
            else:
                s += str(d.age) + " years ago, "
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
        s = self.name + " " + self.title + "\n\n"
        s += "Society type: " + self.society + "\n\n"
        if self.leaderCount == 1:
            s += "Leader: " + self.leader.nameFull() + "\n\n"
        else:
            s += "Leading body: " + self.leader.nameFull() + "\n\n"
        s += "Population: " + str(self.populationCount()) + "\n\n"
        s += "Capital: " + self.origin.city.name + "\n\n"
        return s

# This doesn't strictly represent a "population" as such, 
# just any sort of significant/named entity in the universe.
class Population:
    def __init__(self,c,n=None,t="",a=None,p=1,kind="person"):
        self.parents = []
        self.kids = []
        self.culture = c
        self.number = p
        if n == None or n in self.culture.populations.keys():
            self.name = (self.culture.language.genName(),self.culture.language.genName())
            if self.number > 1:
                self.name = (self.culture.language.genName(),"")
        else:
            self.name = (n,"")
        if a == None:
            self.age = 0
            for i in range(self.number):
                self.age += random.randint(20,40)
            self.age = math.floor(self.age/self.number)
        else:
            self.age = a
        self.title = t + " "
        self.fullName = self.nameFull()
        self.culture.populations[self.name] = self
        self.pronouns = ["it","he","she"]
        self.possessive = ["its","his","her"]
        self.gender = random.choice([0,1,2])
        self.description = ""
        self.kind = kind
        self.associations = []
    def associate(self):
        k = random.choice(list(self.culture.myMap.spheres))
        self.associations.append(k)
        self.descrip()
    def descrip(self):
        s = self.pronouns[self.gender].capitalize() + " is"
        if self.age < self.culture.mythAge:
            s += " a " + str(self.age) + "-year-old "
        else:
            s += " an ageless "
        if self.kind == "deity":
            s += "deity of "
            s += synonym(self.associations[0],seedNum(self.name[0]))
            if self.associations[0] == "mountain":
                s += "s"
        if self.kind == "person":
            s += "person from the "
            s += self.culture.name + " culture"
        if self.kind == "location":
            s += "location of the "
            s += self.culture.name + " culture"
        s += ".\n"
        if len(self.parents) == 0:
            s += " "
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
            s += " "
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
        self.description = s
        return s
    def addPop(self,p):
        self.age = 0
        for g in range(self.number):
            self.age += self.age
        for f in range(p.number):
            self.age += p.age
        self.number = self.number+p.number
        self.age = self.age/self.number
        p.number = 0
    def agePop(self,scl):
        self.age += scl
        if self.age > self.culture.oldAge:
            self.number = round(self.number*random.random())
    def justName(self):
        s = self.name[0]
        if self.name[1] != "":
            s += " " + self.name[1]
        return s
    def nameFull(self):
        if self.title != "":
            s = self.title + self.name[0]
        else:
            s = self.kind + " " + self.name[0]
        if self.name[1] != "":
            s += " " + self.name[1]
        return s
    def popNotes(self):
        s = "The "
        s += self.nameFull()
        s += " is a " + self.kind
        s += " of the society of " + self.culture.name + ".\n\n"
        s += self.descrip()
        return s

class Flag:
    def __init__(self,c):
        self.culture = c
        self.xDim = 384
        self.yDim = 192
        self.colors = []
        self.newColor()
        self.genFlag()
    def newColor(self):
        h = random.randint(0,255)
        s = random.randint(0,255)
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
    def genFlag(self):
        img = Image.new('HSV',(self.xDim,self.yDim),self.colors[0])
        drawer = ImageDraw.Draw(img)
        numElements = random.randint(1,4)
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
        self.filename = "flag_" + self.culture.name + ".gif"
        img = img.convert('RGB')
        img.save(self.filename,"GIF")

class Language:
    def __init__(self,c):
        self.culture = c
        self.characters()
        self.lengthPref = random.choice([3,5,9])
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
    def genName(self):
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
            if k == 0 and random.random() < 0.35:
                c = self.preferredStart
            n += c
            lastchar = c
        return n.capitalize()

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
        self.date = 0
        self.seasonStrength = 0.05
        self.setNorth()
        self.biomeColors()
        self.displayNo = None
        self.infoGui = None
        self.viewmode = 0
        self.timeScale = 1
        self.age = random.randint(1000,100000)
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
        self.north = Node(nx,ny)
    def nodeLat(self,n):
        return "Latitude: " + str(round(n.y/self.distScale))
    def nodeLong(self,n):
        return "Longitude: " + str(round(n.x/self.distScale))
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
            reg = "Inside " + n.resourceRegion.rootCity.name + " outskirts" + "\n"
            reg += "Total food available: " + str(math.floor(n.resourceRegion.resources[0]*self.rscScale)) + " t/year \n"
            reg += "Total industrial resources available: " + str(math.floor(n.resourceRegion.resources[1]*self.rscScale)) + " t/year \n"
        return reg
    def infoScales(self):
        self.distScale = 12
        self.eScale = random.randint(2000,3000)
        self.tempScale = 105
        self.rainfallScale = 7628
        self.fertScale = 100
        self.metalScale = 140000
        self.vegScale = 10295
        self.wildlifeScale = 500
        self.rscScale = 5
    def nodeInfo(self,n):
        self.divWidth = 64
        info = ""
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
        search = Node(xx,yy)
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
        if root.landmass != None or root.elevation < self.sealevel:
            return -1
        root.landmass = Landmass(root,self.sealevel)
        self.landmasses.append(root.landmass)
    def buildAllLand(self):
        print("Building landmasses...")
        for p in self.atlas:
            self.buildLandmass(p)
    def buildBodyWater(self,root):
        if root.landmass != None or root.bodyWater != None:
            return -1
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
        hillCenter = Node(xx,yy)
        for p in self.atlas:
            dist = p.dist(hillCenter)
            if dist <= radius:
                multiplier = distMod(dist,radius)
                if p.elevation > 0.75:
                    multiplier = multiplier*((1-p.elevation)**2)
                p.elevation = p.elevation+(maximum*multiplier)
    def addHill(self,xx,yy,maximum=0.25,radius=128):
        hillCenter = Node(xx,yy)
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
            self.addSineHill(self.xDim/2,self.yDim/2,0.4,radius=self.xDim/1.5)
            if shape == "volcanic":
                self.smooth(4)
                self.elevationAdd(-0.1)
                self.addSineHill(self.xDim/2,self.yDim/2,0.25,radius=self.xDim*1.5)
                self.smooth(4)
                self.addHill(self.xDim/2,self.yDim/2,0.45,radius=self.xDim/11)
                self.addSineHill(self.xDim/2,self.yDim/2,0.-0.1,radius=self.xDim/11)
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
        if shape == "archipelago":
            self.addHills(16,0.25)
        if shape == "plain":
            self.addSineHill(self.xDim/2,self.yDim/2,0.4,radius=self.xDim*5)
        if shape != "volcanic":
            self.addMountains()
            self.addHills()
            self.erode(3)
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
                            n.elevation = n.elevation*0.98
    def nodeSlopes(self):
        for p in self.atlas:
            p.getSlope()
    def rainfall(self):
        for p in self.atlas:
            rainfall = 0.4*random.uniform(0.95,1.05)
            rainfall = ((rainfall*(0.45-((p.temp*0.5)+(p.elevation*1.5))))+rainfall)/2
            if p.river != None:
                rainfall = rainfall*1.5
            for l in p.neighbors:
                if l.elevation < self.sealevel:
                    rainfall *= 1.25
            p.rainfall = clamp(rainfall,0,1)
        for p in self.atlas:
            p.rainfall = sum([n.rainfall for n in p.neighbors])/len(p.neighbors)
    def temperature(self):
        for p in self.atlas:
            s = self.seasonStrength
            seasonMod = clamp(1+(math.sin((math.pi*self.date)/2)*s),1-s,1+s)
            latTemp = 0.375*(p.dist(self.north)/self.xDim)
            elevationTemp = 0.575*(1-p.elevation)
            p.temp = clamp(seasonMod*(elevationTemp+latTemp),0,1)
    def soilProperties(self):
        for p in self.atlas:
            p.metallicity = 0.5*p.elevation*random.uniform(0.666,1.5)
            if p.river != None:
                p.metallicity *= 1.25
            p.metallicity = clamp(p.metallicity,0,1)
            fertilityBase = abs(p.elevation-(self.sealevel*1.2))*random.uniform(0.8,1.25)
            if fertilityBase == 0:
                p.fertility = 1
            else:
                p.fertility = 0.07/fertilityBase
            if p.river != None:
                p.fertility *= 2
            p.fertility = clamp(p.fertility,0,1)
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
            slope = clamp((p.realSlope()*(12000)),-32,32)
            shade = math.floor((-16)+p.biomeColor[2]+slope+(((p.elevation+1)**3)*16))
            p.biomeColor = (p.biomeColor[0],p.biomeColor[1],shade)
    def setWildlife(self):
        for p in self.atlas:
            p.herbivores = clamp((p.vegetation*random.uniform(0.8,1.25))**2,0,1)
            p.carnivores = clamp(((p.herbivores*1.5)**2),0,1)
    def biomeColors(self):
        bColors = {}
        bColors["desert"] = (16,64,142)
        bColors["savanna"] = (48,64,136)
        bColors["tundra"] = (128,48,142)
        bColors["shrubland"] = (72,96,128)
        bColors["boreal forest"] = (110,86,108)
        bColors["forest"] = (96,128,128)
        bColors["tropical forest"] = (78,162,96)
        bColors["frost"] = (134,32,206)
        bColors["mountain"] = (134,0,136)
        bColors["water"] = (142,64,64)
        self.biomeColors = bColors
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
                   "sailors":0.8,
                   "astrology":0.15,
                   "freedom":0.35,
                   "warriors":0.05}
        self.values["food"] = {"agriculture":-0.1,
                   "greed":-0.25,
                   "materialists":0.55,
                   "collectivists":0.25,
                   "simplicity":0.45,
                   "worshippers":-0.2,
                   "superstition":-0.2,
                   "traders":0.2,
                   "warriors":0.2,
                   "builders":0.25}
        self.values["darkness"] = {"travelers":0.2,
                   "collectivists":-0.4,
                   "superstition":0.7,
                   "greed":0.55,
                   "astrology":0.25,
                   "materialists":-0.15,
                   "shamans":0.15,
                   "freedom":-0.2,
                   "warriors":0.5,
                   "worshippers":0.4}
        self.values["movement"] = {"travelers":0.8,
                   "sailors":0.2,
                   "traders":0.55,
                   "astrology":0.15,
                   "builders":-0.4,
                   "simplicity":0.4,
                   "materialists":-0.15,
                   "freedom":0.85,
                   "naturalists":0.2,
                   "collectivists":-0.2,
                   "warriors":0.2,
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
        self.values["nature"] = {"naturalists":0.35,
                   "shamans":0.2,
                   "agriculture":0.25,
                   "freedom":0.25,
                   "travelers":0.15,
                   "simplicity":0.5,
                   "collectivists":0.2,
                   "superstition":0.3,
                   "astrology":0.1,
                   "metallurgists":0.4,
                   "warriors":0.05,
                   "worshippers":0.15}
        self.values["growth"] = {"shamans":0.1,
                   "agriculture":0.55,
                   "naturalists":0.45,
                   "metallurgists":-0.2,
                   "freedom":0.15,
                   "astrology":-0.3,
                   "collectivists":0.4,
                   "materialists":-0.1,
                   "warriors":0.05,
                   "builders":0.15,
                   "greed":0.1,
                   "worshippers":0.1}
        self.values["sky"] = {"travelers":0.55,
                   "craftsmen":0.3,
                   "traders":0.25,
                   "superstition":0.5,
                   "metallurgists":0.1,
                   "worshippers":0.55,
                   "freedom":0.55,
                   "shamans":0.1,
                   "astrology":0.7,
                   "simplicity":0.4,
                   "builders":0.1}
        self.values["earth"] = {"metallurgists":1.2,
                   "craftsmen":0.9,
                   "traders":0.45,
                   "materialists":0.75,
                   "agriculture":0.25,
                   "collectivists":0.4,
                   "builders":0.6,
                   "freedom":-0.25,
                   "worshippers":-0.15,
                   "greed":0.55,
                   "superstition":-0.25,
                   "warriors":0.05,
                   "astrology":-0.2}
        self.values["fields"] = {"agriculture":0.75,
                   "builders":0.3,
                   "materialists":0.5,
                   "naturalists":0.35,
                   "superstition":-0.3,
                   "simplicity":-0.3,
                   "collectivists":0.25,
                   "warriors":0.05}
        self.values["sunlight"] = {"worshippers":0.9,
                   "astrology":0.6,
                   "naturalists":0.15,
                   "travelers":0.3,
                   "simplicity":0.75,
                   "freedom":0.6,
                   "materialists":-0.2,
                   "warriors":0.3,
                   "builders":-0.1}
        self.values["ice"] = {"superstition":0.2,
                   "simplicity":-0.4,
                   "freedom":-0.4,
                   "travelers":-0.45,
                   "materialists":0.45,
                   "sailors":0.1,
                   "shamans":0.45,
                   "greed":0.55,
                   "metallurgists":0.3,
                   "warriors":0.4,
                   "agriculture":-0.2,
                   "collectivists":0.15}
        self.values["fear"] = {"superstition":0.75,
                   "worshippers":0.5,
                   "shamans":0.75,
                   "freedom":-0.5,
                   "collectivists":0.35,
                   "simplicity":-0.3,
                   "builders":0.3,
                   "greed":0.7,
                   "materialists":-0.25,
                   "warriors":0.45}
        self.values["death"] = {"freedom":-0.1,
                   "collectivists":-0.1,
                   "warriors":0.55,
                   "greed":0.15,
                   "travelers":0.15,
                   "builders":-0.1,
                   "simplicity":-0.1,
                   "materialists":0.2}
        self.values["water"] = {"sailors":2,
                   "simplicity":0.3,
                   "freedom":0.75,
                   "travelers":0.6,
                   "builders":0.15,
                   "craftsmen":0.1,
                   "astrology":0.45,
                   "traders":0.75,
                   "warriors":0.1}
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
        self.influences["elevation"] = {"sky":0.3,
                       "sunlight":0.15,
                       "fields":0.1,
                       "earth":0.1,
                       "ice":0.1,
                       "fear":0.05,
                       "movement":0.4,
                       "death":0.1}
        self.influences["hills"] = {"sky":0.1,
                       "sunlight":0.1,
                       "earth":0.1,
                       "fear":0.05,
                       "movement":0.2,
                       "death":0.1}
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
        self.influences["rainfall"] = {"water":0.5,
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
        self.influences["rivers"] = {"water":2.5,
                       "swimming":1.0,
                       "nature":0.1,
                       "food":0.35,
                       "growth":0.1,
                       "movement":0.3}
        self.influences["water"] = {"water":4.0,
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
                       "fear":0.2,
                       "darkness":0.2,
                       "fields":0.15,
                       "water":0.1,
                       "death":0.4}
    def godSpheres(self):
        s0 = list(self.influences.keys())
        s1 = list(self.values.keys())
        s2 = list(self.valuesOutputs.keys())
        self.spheres = s0+s1+s2
    def nearestCityDist(self,xx,yy):
        if len(self.cities) < 1:
            return self.xDim
        else:
            c = self.nearestCity(xx,yy).node
            d = Node(xx,yy).dist(c)
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
        while (cityNode.biome == "water" or cityNode.city != None or
               cityNode.x < 32 or cityNode.x > self.xDim-32 or cityNode.y < 32 
               or cityNode.y > self.yDim-32 or self.nearestCityDist(cityNode.x,cityNode.y) < 32):
            cityNode = random.choice(self.atlas)
        newCity = City(cityNode,pop=random.randint(12,136),m=self)
    def scatterCities(self,n):
        for k in self.atlas:
            k.defaultRoads()
        for i in range(n):
            self.randomCity()
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
        clickedNode = self.nearestNode(event.x,event.y)
        self.displayString.set(self.nodeInfo(clickedNode))
        self.displayNo = clickedNode
        cityNode = self.nearestCity(event.x,event.y)
        if cityNode.node.dist(Node(event.x,event.y)) < 8:
            self.displayString.set(self.nodeInfo(cityNode.node))
            self.displayNo = cityNode.node
    def redraw(self):
        visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
        graphDraw = ImageDraw.Draw(visualAtlas)
        if self.viewmode == 0:
            for tri in self.triangles:
                tri.drawReal(graphDraw,self.sealevel)
            for n in self.atlas:
                n.drawReal(graphDraw,self.sealevel)
            for l in self.landmasses:
                for r in l.rivers:
                    r.drawRiver(graphDraw,self.xDim)
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
        self.setBiomes()
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
        for c in self.cities:
            c.updateDemog()
    def updatePops(self):
        for c in self.cultures:
            c.updatePops()
    def updateEvents(self):
        for e in self.events:
            e.ageEvent()
    def nextTurn(self):
        self.updateTiming()
        self.updateResources()
        self.updateDemogs()
        self.updatePops()
        self.updateEvents()
        self.updateTerritory()
        self.redraw()
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
        b1 = Button(self.infoGui,text="Society Info",command=self.cultureInfo)
        b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
        c1 = "medium aquamarine"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        if p.kind == "deity":
            b1 = Button(self.infoGui,text="Mythology Info",command=self.mythologyInfo)
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "light goldenrod"
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
            c1 = "SteelBlue4"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def cultureInfo(self):
        if self.displayNo == None:
            return -1
        if self.displayNo.culture == None:
            return -1
        self.displayCulture = self.displayNo.culture
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        photo = Image.open(self.displayCulture.flag.filename)
        self.flagImg = ImageTk.PhotoImage(photo)
        self.flagLbl = Label(self.infoGui,image=self.flagImg)
        self.flagLbl.config(borderwidth=32)
        self.flagLbl.photo = photo
        self.flagLbl.pack()
        self.cultureString = StringVar()
        self.cultureString.set(self.displayCulture.cultureNotes())
        cdsc = Label(self.infoGui,textvariable=self.cultureString)
        cdsc.pack()
        b1 = Button(self.infoGui,text="Mythology Info",command=self.mythologyInfo)
        b1.pack(anchor=S,side=RIGHT)
        c1 = "light goldenrod"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def mythologyInfo(self):
        if self.displayCulture == None:
            return -1
        if self.infoGui != None:
            self.infoGui.destroy()
        self.infoGui = Toplevel()
        self.mythString = StringVar()
        self.mythString.set(self.displayCulture.mythNotes())
        mdsc = Label(self.infoGui,textvariable=self.mythString)
        mdsc.config(justify=LEFT)
        mdsc.pack(anchor=W,side=RIGHT)
        for i in range(len(self.displayCulture.deities)):
            s = " "+self.displayCulture.deities[i].justName()+" "
            b1 = Button(self.infoGui,text=s)
            b1.configure(command = lambda self=self, d = self.displayCulture.deities[i]: self.popInfo(d))
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "SteelBlue4"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def cityInfo(self):
        if self.displayNo == None:
            return -1
        if self.displayNo.city == None:
            return -1
        self.displayCity = self.displayNo.city
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
        self.cityString.set(self.displayCity.cityNotes())
        cdsc = Label(self.infoGui,textvariable=self.cityString)
        cdsc.pack()
        self.displayCulture = self.displayNo.culture
        b1 = Button(self.infoGui,text="Society Info",command=self.cultureInfo)
        b1.pack(anchor=S,side=RIGHT)
        c1 = "medium aquamarine"
        b1.config(bg=c1,activebackground=c1,activeforeground=c1)
    def changeView(self):
        if self.viewmode == 1:
            self.viewmode = 0
        else:
            self.viewmode += 1
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
            self.mapname = "map_" + self.cultures[0].language.genName() + ".gif"
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
            b1 = Button(gui,text="Society Info",command=self.cultureInfo)
            b1.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "medium aquamarine"
            b1.config(bg=c1,activebackground=c1,activeforeground=c1)
            b3 = Button(gui,text="Settlement Info",command=self.cityInfo)
            b3.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "aquamarine"
            b3.config(bg=c1,activebackground=c1,activeforeground=c1)
            b2 = Button(gui,text="Change Mode",command=self.changeView)
            b2.pack(anchor=S,side=TOP,expand=YES,fill=BOTH)
            c1 = "sandy brown"
            b2.config(bg=c1,activebackground=c1,activeforeground=c1)
            gui.mainloop()

#----------------------------------------------------------------------#            
# Let's generate a map

numNodes = 2**14
mapDimX = 960
mapDimY = 960
atlas = [Node(-1,-1),Node(mapDimX+1,-1),Node(mapDimY+1,mapDimY+1),Node(-1,mapDimY+1)]
world = Map(atlas,numNodes,mapDimX,mapDimY)

print("Generating points...")
for x in range(numNodes-4):
    nodeX = random.random()*mapDimX
    nodeY = random.random()*mapDimY
    newNode = Node(nodeX,nodeY)
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
relaxLloyd(npFloatAtlas,1)
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
    newTri = Triangle(atlas[triVertsIndices[0]],atlas[triVertsIndices[1]],atlas[triVertsIndices[2]])
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
world.buildAllLand()
world.buildAllWater()
world.addMajorRiver(12)
world.addMinorRiver(12)
world.cullStreams()
world.soilProperties()
print("Defining biomes...")
world.setBiomes()
world.buildRegions()
world.setWildlife()
world.influences()
world.values()
world.godSpheres()
world.scatterCities(16)
print("Drawing map...")
root = Tk()
world.drawReal(root)