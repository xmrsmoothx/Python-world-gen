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
    def setVegetation(self):
        tempFitness = 1-abs(0.3-self.temp)
        elevationFitness = 1-abs(0.45-self.elevation)
        fertilityFitness = self.fertility+(self.metallicity*0.25)
        rainFitness = 1-abs(self.rainfall-0.8)
        vegFitness = ((tempFitness+1)*(elevationFitness+0.5)*(fertilityFitness+0.5)*(rainFitness+2))
        self.vegetation = clamp(vegFitness/16,0,1)
    def setBiome(self,sl):
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
    def claim(self,n):
        n.culture = self.culture
        n.allegiance = self.allegiance+1
        if self.culture.name not in n.region.culturalNames:
            n.region.culturalNames[self.culture.name] = self.culture.language.genName()
        if n.landmass != None:
            if self.culture.name not in n.landmass.culturalNames:
                n.landmass.culturalNames[self.culture.name] = self.culture.language.genName()
        if n.river != None:
            if self.culture.name not in n.river.culturalNames:
                n.river.culturalNames[self.culture.name] = self.culture.language.genName()
    def updateAllegiance(self):
        if self.culture != None:
            chance = 1/clamp(self.allegiance,0.00001,512)
            for n in self.neighbors:
                roll = random.random()
                if roll <= chance:
                    self.claim(n)
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
        if s <= 3:
            t = "isle"
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
        inf = "slope"
        self.setOutput(inf,strength)
        strength = self.node.vegetation
        inf = "vegetation"
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
    def cType(self,p):
        if p <= 20:
            c = random.choice(["bivouac","camp","camp","encampment","campsite"])
        elif p <= 100:
            c = random.choice(["village","village","hamlet"])
        elif p <= 500:
            c = random.choice(["township","settlement"])
        elif p <= 2000:
            c = "town"
        elif p <= 10000:
            c = "city"
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
        col = (0,0,255)
        out = (0,0,0)
        if self.population <= 20:
            self.drawTent(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 100:
            self.drawHut(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 500:
            self.drawVillage(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 2000:
            self.drawTown(drawer,self.node.x,self.node.y,col,out)
        elif self.population <= 10000:
            self.drawCity(drawer,self.node.x,self.node.y,col,out)

class ResourceRegion:
    def __init__(self,c,m):
        self.rootCity = c
        self.myMap = m
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
        self.resources[0] = (rawPlant + (rawAnimal*0.7))*self.totalPop
        self.resources[1] = (rawMetal + (rawAnimal*0.3) + (rawPlant*0.4))*self.totalPop
        if "craftsmen" in m:
            self.resources[1] *= 1.1
        if "builders" in m:
            self.resources[1] *= 1.1
        if "collectivists" in m:
            self.resources[0] *= 1.05
    def updateReg(self):
        for p in self.nodes:
            if p.city != None:
                p.resourceDist = math.log2(p.city.population)
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
        self.title = self.setTitle()
    def setSociety(self):
        m = self.value.mainValues
        if "greed" in m and "worshippers" in m and "builders" in m:
            return "Empire"
        if "warriors" in m and "greed" in m and "builders" in m and ("travelers" in m or "sailors" in m):
            return "Imperium"
        if "warriors" in m and "collectivists" in m and "worshippers" in m:
            return "Hegemony"
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
        if "worshippers" in m and "warriors" in m and ("superstitious" in m or "collectivists" in m):
            return "Religious zealots"
        if "shamans" in m and ("naturalists" in m or "astrologists" in m) and "superstitious" in m:
            return "Shamanistic tribe"
        if "shamans" in m and "warriors" in m and ("astrologists" in m or "superstitious" in m or "worshippers" in m):
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
        if ("craftsmen" in m or "metallurgy" in m or "builders" in m) and "simplicity" in m and ("superstitious" in m or "astrologists" in m):
            return "Traditionalist artisans"
        if "greed" in m and "sailors" in m and "warriors" in m:
            return "Pirates"
        if ("travelers" in m or "sailors" in m) and "greed" in m and "warriors" in m:
            return "Raiders"
        if ("travelers" in m or "sailors" in m) and "greed" in m and "simplicity" in m:
            return "Scavengers"
        if "travelers" in m and "simplicity" in m and "freedom" in m:
            return "Hunter-gatherer tribe"
        if "astrologists" in m and "superstitious" in m and "worshippers" in m:
            return "Religious sovereignty"
        if "collectivists" in m and "agriculture" in m and "naturalists" in m:
            return "Agriculturalists"
        if "travelers" in m and "simplicity" in m and ("metallurgists" in m or "craftsmen" in m or "builders" in m):
            return "Nomadic artisans"
        if "travelers" in m and "simplicity" and ("naturalists" in m or "superstitious" in m or "astrologists" in m or "shamans" in m):
            return "Nomadic tribe"
        if "materialists" in m and "metallurgists" in m and ("craftsmen" in m or "builders" in m):
            return "Blacksmiths"
        if "warriors" in m and "collectivists" in m and ("travelers" in m or "sailors" in m):
            return "Revolutionary commune"
        if "builders" in m and "metallurgists" in m and "craftsmen" in m and "agriculture" in m:
            return "Syndicalists"
        if "warriors" in m and "builders" in m and ("worshippers" in m or "superstitious" in m) and "freedom" not in m and "traders" not in m:
            return "Nationalists"
        return "Tribe"
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
            self.society == "Traders" or self.society == "Independent merchants" or self.society == "Mercantile folk"):
            return "Caravans"
        if (self.society == "Blacksmiths" or self.society == "Traditionalist artisans" or
            self.society == "Naturalist artisans" or self.society == "Cooperative artisans"):
            return "Artisans"
        if (self.society == "Socialists" or self.society == "Syndicalists" or self.society == "Revolutionary commune"):
            return random.choice(["People's Union","Union","Collective"])
        if (self.society == "Shamanistic warriors" or self.society == "Shamanistic tribe"):
            return "Mystics"
        if self.society == "Pirates" or self.society == "Raiders":
            return "Brigands"
        return "People"
    def shortName(self):
        name = ""
        name += self.name + " " + self.title
        return name
    def information(self):
        info = ""
        info += self.name +" "+ self.title + "\n"
        info += "("+self.society+")" + "\n"
        return info

class Language:
    def __init__(self,c):
        self.culture = c
        self.characters()
        self.lengthPref = random.choice([3,5,9])
        self.name = self.genName()
    def characters(self):
        c = ['b','c','d','f','g','h','j','k','l','m','n','p','q','r','s','t','v','w','x','y','z']
        v = ['a','e','i','o','u']
        self.langConsonants = []
        self.langVowels = []
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
        self.resourceRegions = []
        self.sealevel = 0.4
        self.setNorth()
        self.biomeColors()
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
                q += n.region.biome
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
            reg += "Area food output: " + str(math.floor(n.resourceRegion.resources[0]*self.resourceScale)) + " t/year \n"
            reg += "Area industrial output: " + str(math.floor(n.resourceRegion.resources[1]*self.resourceScale)) + " t/year \n"
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
    def nodeInfo(self,n):
        self.divWidth = 48
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
        n = self.cities[0]
        minDist = 1000000
        search = Node(xx,yy)
        for p in self.cities:
            dist = p.node.dist(search)
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
                self.addSineHill(xx,yy,0.3,radius=self.xDim*2)
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
            rainfall = ((rainfall*(0.4-(p.temp+p.elevation)))+rainfall)/2
            q = p
            distance = self.n/256
            diff = 0
            for i in range(math.floor(distance)):
                q = p.westMostNeighbor()
                diff += q.elevation-p.elevation
            diff = diff/distance
            if diff > 0.2 and 0:
                rainfall = clamp(rainfall*(0.5-diff),0,1)
            if p.river != None:
                rainfall = rainfall*1.5
            for l in p.neighbors:
                if l.elevation < self.sealevel:
                    rainfall *= 1.25
            p.rainfall = clamp(rainfall,0,1)
    def temperature(self):
        for p in self.atlas:
            latTemp = 0.375*(p.dist(self.north)/self.xDim)
            elevationTemp = 0.575*(1-p.elevation)
            p.temp = clamp(elevationTemp+latTemp,0,1)
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
            slope = clamp(p.realSlope()*(2**16),-48,48)
            shade = math.floor((-16)+p.biomeColor[2]+slope+(((p.elevation+1)**2)*16))
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
        bColors["boreal"] = (110,86,108)
        bColors["forest"] = (96,128,128)
        bColors["tropical"] = (78,162,96)
        bColors["frozen"] = (134,32,206)
        bColors["mountain"] = (134,0,136)
        bColors["water"] = (142,64,64)
        self.biomeColors = bColors
    def values(self):
        self.valuesOutputs = {"travelers":0,
                              "craftsmen":0,
                              "traders":0,
                              "superstitious":0,
                              "metallurgists":0,
                              "worshippers":0,
                              "freedom":0,
                              "shamans":0,
                              "astrologists":0,
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
                   "astrologists":0.15,
                   "freedom":0.35,
                   "warriors":0.05}
        self.values["food"] = {"agriculture":-0.1,
                   "greed":-0.25,
                   "materialists":0.55,
                   "collectivists":0.25,
                   "simplicity":0.45,
                   "worshippers":-0.2,
                   "superstitious":-0.2,
                   "traders":0.2,
                   "warriors":0.2,
                   "builders":0.25}
        self.values["darkness"] = {"travelers":0.2,
                   "collectivists":-0.4,
                   "superstitious":0.7,
                   "greed":0.4,
                   "astrologists":0.25,
                   "materialists":-0.15,
                   "shamans":0.15,
                   "freedom":-0.2,
                   "warriors":0.5}
        self.values["movement"] = {"travelers":0.8,
                   "sailors":0.2,
                   "traders":0.55,
                   "astrologists":0.15,
                   "builders":-0.4,
                   "simplicity":0.4,
                   "materialists":-0.15,
                   "freedom":0.85,
                   "naturalists":0.2,
                   "collectivists":-0.2,
                   "warriors":0.2}
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
                   "collectivists":0.15,
                   "superstitious":0.3,
                   "astrologists":0.1,
                   "metallurgists":0.4,
                   "warriors":0.05}
        self.values["growth"] = {"shamans":0.1,
                   "agriculture":0.55,
                   "naturalists":0.45,
                   "metallurgists":-0.2,
                   "freedom":0.15,
                   "astrologists":-0.3,
                   "collectivists":0.35,
                   "materialists":-0.1,
                   "warriors":0.05,
                   "builders":0.15}
        self.values["sky"] = {"travelers":0.55,
                   "craftsmen":0.3,
                   "traders":0.25,
                   "superstitious":0.5,
                   "metallurgists":0.1,
                   "worshippers":0.5,
                   "freedom":0.55,
                   "shamans":0.1,
                   "astrologists":0.7,
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
                   "greed":0.6,
                   "superstitious":-0.25,
                   "warriors":0.05,
                   "astrologists":-0.2}
        self.values["fields"] = {"agriculture":0.75,
                   "builders":0.3,
                   "materialists":0.5,
                   "naturalists":0.35,
                   "superstitious":-0.3,
                   "simplicity":-0.3,
                   "collectivists":0.25,
                   "warriors":0.05}
        self.values["sunlight"] = {"worshippers":0.9,
                   "astrologists":0.6,
                   "naturalists":0.15,
                   "travelers":0.25,
                   "simplicity":0.75,
                   "freedom":0.6,
                   "materialists":-0.2,
                   "warriors":0.3,
                   "builders":-0.1}
        self.values["ice"] = {"superstitious":0.3,
                   "simplicity":-0.4,
                   "freedom":-0.4,
                   "travelers":-0.45,
                   "materialists":0.4,
                   "sailors":0.1,
                   "shamans":0.45,
                   "greed":0.7,
                   "metallurgists":0.2,
                   "warriors":0.4,
                   "agriculture":-0.2}
        self.values["fear"] = {"superstitious":0.85,
                   "worshippers":0.3,
                   "shamans":0.8,
                   "freedom":-0.55,
                   "collectivists":0.4,
                   "simplicity":-0.3,
                   "builders":0.25,
                   "greed":0.8,
                   "materialists":-0.25,
                   "warriors":0.65}
        self.values["water"] = {"sailors":2,
                   "simplicity":0.3,
                   "freedom":0.75,
                   "travelers":0.6,
                   "builders":0.15,
                   "craftsmen":0.1,
                   "astrologists":0.45,
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
                       "fear":0.05}
        self.influences["slope"] = {"sky":0.1,
                       "sunlight":0.1,
                       "earth":0.1,
                       "fear":0.05}
        self.influences["temperature"] = {"sunlight":0.2,
                       "darkness":-0.1,
                       "sky":0.1,
                       "movement":0.2,
                       "fear":0.15,
                       "water":0.05,
                       "food":0.1}
        self.influences["vegetation"] = {"nature":0.15,
                       "food":0.15,
                       "darkness":0.05,
                       "plantlife":0.2,
                       "sunlight":-0.1,
                       "fields":-0.1}
        self.influences["rainfall"] = {"water":0.5,
                       "growth":0.2,
                       "nature":0.1,
                       "sky":0.25,
                       "movement":0.1,
                       "darkness":0.1,
                       "plantlife":0.1,
                       "swimming":0.15}
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
                       "nature":0.1}
        self.influences["carnivores"] = {"darkness":0.35,
                       "food":0.25,
                       "fear":0.4,
                       "nature":0.15}
        self.influences["herbivores"] = {"nature":0.35,
                       "food":0.65,
                       "growth":0.3,
                       "earth":0.1,
                       "plantlife":-0.1}
        self.influences["river"] = {"water":2.5,
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
                       "fields":0.1}
        self.influences["desert"] = {"earth":0.4,
                       "sky":0.3,
                       "sunlight":0.7,
                       "food":-0.3,
                       "fields":0.35}
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
                       "ice":0.5}
        self.influences["mountain"] = {"food":-0.3,
                       "fields":0.1,
                       "sky":0.7,
                       "sunlight":0.35,
                       "earth":0.45,
                       "fear":0.2,
                       "nature":0.15,
                       "ice":0.5}
        self.influences["tropical"] = {"food":0.35,
                       "fear":0.4,
                       "darkness":0.35,
                       "fields":-0.25,
                       "nature":0.5,
                       "plantlife":0.6,
                       "growth":0.45,
                       "earth":0.1}
        self.influences["boreal"] = {"food":0.1,
                       "nature":0.25,
                       "growth":0.1,
                       "fields":0.2,
                       "darkness":0.2,
                       "plantlife":0.25,
                       "earth":0.15,
                       "ice":0.15}
        self.influences["frozen"] = {"ice":0.75,
                       "earth":0.2,
                       "food":-0.3,
                       "plantlife":-0.2,
                       "fear":0.2,
                       "darkness":0.2,
                       "fields":0.15,
                       "water":0.1}
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
        cityNode = self.atlas[math.floor(random.random()*len(self.atlas))]
        while (cityNode.biome == "water" or cityNode.city != None or
               cityNode.x < 0 or cityNode.x > self.xDim or cityNode.y < 0 or cityNode.y > self.yDim):
            if len(self.cities) >= 1:
                cityNode = self.atlas[math.floor(random.random()*len(self.atlas))]
                while cityNode.dist(self.nearestCity(cityNode.x,cityNode.y).node) < 32:
                    cityNode = self.atlas[math.floor(random.random()*len(self.atlas))]
            else:
                cityNode = self.atlas[math.floor(random.random()*len(self.atlas))]
        newCity = City(cityNode,pop=random.randint(12,136),m=self)
    def scatterCities(self,n):
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
        cityNode = self.nearestCity(event.x,event.y)
        if cityNode.node.dist(Node(event.x,event.y)) < 8:
            self.displayString.set(self.nodeInfo(cityNode.node))
    def redraw(self):
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
        visualAtlas = visualAtlas.convert("RGB")
        visualAtlas.save("map00.gif","GIF")
        photo = Image.open("map00.gif")
        self.img = ImageTk.PhotoImage(photo)
        self.lbl.configure(image = self.img)
        self.lbl.image = self.img
    def updateResources(self):
        for r in self.resourceRegions:
            r.updateReg()
    def updateTerritory(self):
        for c in self.cities:
            c.node.culture = c.culture
            c.node.allegiance = 1/c.population
        for p in self.atlas:
            p.updateAllegiance()
    def nextTurn(self):
        self.updateResources()
        self.updateTerritory()
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
            desc.pack(side=RIGHT,fill=Y)
            visualAtlas = visualAtlas.convert("RGB")
            visualAtlas.save("map00.gif","GIF")
            photo = Image.open("map00.gif")
            self.img = ImageTk.PhotoImage(photo)
            self.lbl = Label(gui,image=self.img)
            self.lbl.pack()
            self.lbl.bind("<Button-1>",self.displayNode)
            b = Button(gui,text="Next Turn",command=self.nextTurn)
            b.pack(anchor=E,side=RIGHT)
            self.nextTurn()
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
print("Defining biomes...")
world.setBiomes()
world.buildRegions()
world.setWildlife()
world.influences()
world.values()
world.scatterCities(16)
print("Drawing map...")
root = Tk()
world.drawReal(root)