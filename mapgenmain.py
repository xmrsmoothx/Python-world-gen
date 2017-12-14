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
        tempFitness = 1-abs(self.temp-0.7)
        elevationFitness = clamp((1-self.elevation)-(self.slope/2),0,1)
        fertilityFitness = self.fertility-self.metallicity
        rainFitness = self.rainfall
        vegFitness = (tempFitness*0.2+elevationFitness*0.2+fertilityFitness*0.3+rainFitness*0.25)
        self.vegetation = clamp(vegFitness,0,1)
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
            if self.rainfall < 0.07:
                self.biome = "desert"
            else:
                self.biome = "tropical"
        if self.elevation < sl:
            self.biome = "water"
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
    def removeRiver(self):
        for n in self.nodes:
            self.nodes.remove(n)
            n.river = None
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
            scale = xDim/4
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
    def centermass(self):
        xTotal = sum([p.x for p in self.nodes])
        yTotal = sum([p.y for p in self.nodes])
        xx = xTotal/self.size
        yy = yTotal/self.size
        self.centroid = Node(xx,yy)
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

class influence:
    def __init__(self,myMap,myInfluence):
        self.influenceType = myInfluence
        self.setOutput(myMap)
    def setOutput(self,myMap):
        self.influenceOutput = {}
        for o in myMap.influenceOutputs:
            if myMap.influences[self.influenceType].has_key(o):
                self.influenceOutput[o] = myMap.influences[self.influenceType][o]
            else:
                self.influenceOutput[o] = 0

class Map:
    def __init__(self,aAtlas,numNodes,mapDimX,mapDimY):
        self.atlas = aAtlas
        self.n = numNodes
        self.xDim = mapDimX
        self.yDim = mapDimY
        self.landmasses = []
        self.waterBodies = []
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
        if n.biome != "ocean":
            return "Temperature: " + str(round((n.temp*self.tempScale)-30,1)) + " degrees"
        else:
            return "Temperature: " + str(round((n.temp*self.tempScale*0.3),1)) + " degrees"
    def nodeRain(self,n):
        return "Rainfall: " + str(round(n.rainfall*self.rainfallScale,1)) + "cm/yr"
    def nodeBiome(self,n):
        return n.biome
    def infoScales(self):
        self.distScale = 12
        self.eScale = 2000
        self.tempScale = 105
        self.rainfallScale = 257
    def nodeInfo(self,n):
        info = "          Location Information:          \n"
        info += self.nodeLat(n) + "\n"
        info += self.nodeLong(n) + "\n"
        info += self.nodeElevation(n) + "\n"
        info += self.nodeTemp(n) + "\n"
        info += self.nodeRain(n) + "\n"
        info += self.nodeBiome(n) + "\n"
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
    def addMountains(self,num=5,height=0.25):
        avgRad = self.xDim/3.5
        for i in range(num):
            xx = math.floor(random.random()*self.xDim)
            yy = math.floor(random.random()*self.yDim)
            hillRad = avgRad*random.uniform(0.8,1.25)
            hillHeight = height*random.uniform(0.8,1.25)
            self.addSineHill(xx,yy,hillHeight,hillRad)
    def addHills(self,num=8,height=0.1):
        avgRad = self.xDim/5
        for i in range(num):
            xx = math.floor(random.random()*self.xDim)
            yy = math.floor(random.random()*self.yDim)
            hillRad = avgRad*random.uniform(0.8,1.25)
            hillHeight = height*random.uniform(0.8,1.25)
            self.addSineHill(xx,yy,hillHeight,hillRad)
    def addShape(self,shape):
        if shape == "island":
            self.addSineHill(self.xDim/2,self.yDim/2,0.4,radius=self.xDim/1.5)
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
        if shape == "archipelago":
            self.addHills(16,0.2)
        if shape == "plain":
            self.addSineHill(self.xDim/2,self.yDim/2,0.4,radius=self.xDim*5)
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
            fertilityBase = abs(p.elevation-(self.sealevel*1.15))
            if fertilityBase == 0:
                p.fertility = 1
            else:
                p.fertility = 1/fertilityBase
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
    def influences(self):
        self.influenceOutputs = {"sky":1,
                                 "sunlight":1,
                                 "fields":1,
                                 "earth":1,
                                 "ice":1,
                                 "fear":1,
                                 "water":1,
                                 "growth":1,
                                 "nature":1,
                                 "plantlife":1,
                                 "movement":1,
                                 "darkness":1,
                                 "swimming":1,
                                 "food":1}
        self.influences = {}
        self.influences["elevation"] = {"sky":0.5,
                       "sunlight":0.2,
                       "fields":0.1,
                       "earth":0.1,
                       "ice":0.1,
                       "fear":0.05}
        self.influences["slope"] = {"sky":0.2,
                       "sunlight":0.1,
                       "earth":0.3,
                       "fear":0.1}
        self.influences["rainfall"] = {"water":0.35,
                       "growth":0.2,
                       "nature":0.1,
                       "sky":0.25,
                       "movement":0.1,
                       "darkness":0.1,
                       "plantlife":0.1}
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
                       "food":0.5,
                       "growth":0.3,
                       "earth":0.1}
        self.influences["river"] = {"water":0.5,
                       "swimming":0.2,
                       "nature":0.1,
                       "food":0.35,
                       "growth":0.1,
                       "movement":0.3}
        self.influences["water"] = {"water":0.6,
                       "swimming":0.3,
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
                       "fields":0.15}
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
            img = ImageTk.PhotoImage(photo)
            lbl = Label(gui,image=img)
            lbl.pack()
            lbl.bind("<Button-1>",self.displayNode)
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
world.addShape("archipelago")
world.addMountains()
world.addHills()
world.setSeaLevel(0.4)
world.erode(3)
world.smooth(3)
world.cullDots()
world.clampElevation()
world.buildAllLand()
world.buildAllWater()
world.addMajorRiver(12)
world.addMinorRiver(12)
world.cullStreams()
print("Defining biomes...")
world.setBiomes()
world.setWildlife()
world.influences()
print("Drawing map...")
root = Tk()
world.drawReal(root)