# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 14:52:37 2017

@author: Bri
"""

import numpy as np
import random
import math
from PIL import Image, ImageDraw, ImageTk
from scipy.spatial import Voronoi
import matplotlib.path as mpltPath
from src_tools import *

XDIM = 720
BORDERSCALE = 1/12
STREETRANGE = 6400

def chaoticFunction(x):
    return (23*(math.sin(math.cos(x))*math.sin(x/3)*math.cos(x/5)))+(1*math.sin(x/3))

#  Returns whether it crosses the x axis at or around here
def sinRound(x):
    if abs(chaoticFunction(x)) < 0.6:
        return True
    return False

def xRange(name):
    start = (seedNum(name)*73 % 250)
    end = start+STREETRANGE
    return (start,end)

def yRange(name):
    start = (seedNum(name)*113 % 250)
    end = start+STREETRANGE
    return (start,end)

def xStreets(x,name):
    xStart = xRange(name)[0]
    xReal = xStart+((x*STREETRANGE)/XDIM)
    if sinRound(xReal):
        return True
    return False

def yStreets(y,name):
    yStart = yRange(name)[0]
    yReal = yStart+((y*STREETRANGE)/XDIM)
    if sinRound(yReal):
        return True
    return False

def clamp(x,minimum,maximum):
    if x < minimum:
        return minimum
    elif x > maximum:
        return maximum
    else:
        return x

def exclus(g):
    if g < 0.01 and g > 0:
        return 0.01
    elif g > -0.01 and g < 0:
        return -0.01
    return g

def drawCircle(drawer,x,y,radius,color):
    x1 = x-radius
    x2 = x+radius
    y1 = y-radius
    y2 = y+radius
    drawer.ellipse([(x1,y1),(x2,y2)],color)

def enorm(a,b):
    j = abs(a**2) + abs(b**2)
    return math.sqrt(j)

def colAvg(c0,c1):
    cc = (math.floor((c0[0]+c1[0])/2),math.floor((c0[1]+c1[1])/2),math.floor((c0[2]+c1[2])/2))
    return cc

class StreetNode:
    def __init__(self,coords):
        self.x = coords[0]
        self.y = coords[1]
        self.neighbors = []
        self.type = None
        self.drawColor = (0,0,0)
        self.road = False
        self.node = None
    def dist(self,x,y):
        dx = abs(self.x-x)
        dy = abs(self.y-y)
        dist = math.sqrt((dx**2) + (dy**2))
        return dist
    def nodeDist(self,n):
        x = n.x
        y = n.y
        dx = abs(self.x-x)
        dy = abs(self.y-y)
        dist = math.sqrt((dx**2) + (dy**2))
        return dist
    def midpt(self,n):
        xx = (self.x + n.x)/2
        yy = (self.y + n.y)/2
        return (xx,yy)
    def edgept(self,nbr,seed):
        d = self.dist(nbr.x,nbr.y)
        if d == 0:
            d = 1
        dt = (seed*self.x*self.y*37.13529) % d
        t = dt/d
        xx = ((1-t)*self.x) + (t*nbr.x)
        yy = ((1-t)*self.y) + (t*nbr.y)
        n = StreetNode((xx,yy))
        return n

class Street:
    def __init__(self,n0,n1,w=2):
        self.nodes = [None,None]
        self.nodes[0] = n0
        self.nodes[1] = n1
        n0.neighbors.append(n1)
        n1.neighbors.append(n0)
        self.exists = 1
        self.width = w
    def streetDist(self,x,y):
        d = (self.nodes[0].dist(x,y)+self.nodes[1].dist(x,y))/2
        return abs(d)
    def length(self):
        return self.nodes[0].dist(self.nodes[1].x,self.nodes[1].y)
    def cull(self):
        self.nodes[1].x = self.nodes[0].x
        self.nodes[1].y = self.nodes[0].y
        self.exists = 0
    def drawSelf(self,drawer,col):
        drawCircle(drawer,self.nodes[0].x,self.nodes[0].y,self.width-1,col)
        drawer.line([(self.nodes[0].x,self.nodes[0].y),(self.nodes[1].x,self.nodes[1].y)],fill=col,width=self.width*2)
        drawCircle(drawer,self.nodes[1].x,self.nodes[1].y,self.width-1,col)

class Block:
    def __init__(self,m):
        self.myTown = m
        self.verts = []
        self.neighbors = []
        self.subblocks = []
        self.substreets = []
        self.type = None
        self.node = None
        self.col = (0,0,0)
    def centroid(self):
        if len(self.verts) != 0:
            xx = sum([n.x for n in self.verts])/len(self.verts)
            yy = sum([n.y for n in self.verts])/len(self.verts)
        else:
            xx = 0
            yy = 0
        self.x = xx
        self.y = yy
        self.orderccw()
        return [self.x,self.y]
    def blockDist(self,x,y):
        self.centroid()
        dx = self.x-x
        dy = self.y-y
        dist = math.sqrt((dx*dx) + (dy*dy))
        return dist
    def neighborize(self,j):
        if j not in self.neighbors:
            self.neighbors.append(j)
        if self not in j.neighbors:
            j.neighbors.append(j)
    def sharedNeighbors(self,f):
        s = 0
        for p in self.verts:
            if p in f.verts:
                s += 1
        return s
    def orderccw(self):
        self.verts.sort(key=lambda a: math.atan2(float(self.x-a.x),float(self.y-a.y)))
    def drawRoads(self,drawer):
        dCol = Tools.streetColor
        if len(self.verts) < 3:
            return -1
        vts = [(p.x,p.y) for p in self.verts]
        vts.append(vts[0])
        dCol = Tools.streetColor
        drawer.line(vts,dCol,3)
    def drawBuilding(self,drawer):
        self.drawRoads(drawer)
        dCol = self.col
        if len(self.verts) > 1:
            vts = [(p.x,p.y) for p in self.verts]
            dCol = Tools.streetColor
            bCol = Tools.buildingColor
            drawer.polygon(vts,bCol,dCol)
    def drawSelf(self,drawer):
        dCol = self.col
        if len(self.verts) > 1:
            vts = [(p.x,p.y) for p in self.verts]
            if self.blockDist(self.myTown.x,self.myTown.y) < self.myTown.radius and self.col != Tools.waterColor:
                dCol = Tools.streetColor
                bCol = Tools.buildingColor
                drawer.polygon(vts,bCol,dCol)
            else:
                drawer.polygon(vts,dCol,dCol)

class Town:
    def __init__(self,n,m,nom):
        self.myMap = m
        self.node = n
        self.name = nom
        self.mapName = "./generated/town_"+self.name+".gif"
        self.xDim = XDIM
        self.yDim = XDIM
        #self.landColor = self.myMap.biomeColors[self.node.biome]
        self.landColor = self.node.biomeColor
        self.streetColor = Tools.streetColor
        self.waterColor = Tools.waterColor
        self.x = self.xDim/2
        self.y = self.yDim/2
        self.roadRadius = 3
        if self.node.city != None:
            self.population = self.node.city.population
            self.wateryNeighbors = [u for u in self.node.neighbors if u.watery() == 1]
            self.radius = math.floor(22+(self.population**(17/40)))*(1+(len(self.wateryNeighbors)/10))
            if self.node.river != None:
                self.radius = self.radius+16
        else:
            self.population = 0
            self.radius = 0
        cnt = 1024
        sd = (self.node.x*73)+(self.node.y*37)
        random.seed(sd)
        buffer = 32
        verts = [[random.randint(buffer,self.xDim-buffer),random.randint(buffer,self.yDim-buffer)] for i in range(cnt)]
        verts.append([self.xDim/2,self.yDim/2])
        gridSize = 20
        for p in verts:
            if random.random() < 0.5:
                p[0] = round(p[0]/gridSize)*gridSize
                p[1] = round(p[1]/gridSize)*gridSize
            else:
                if random.random() < 0.5:
                    p[1] = round(p[1]/gridSize)*gridSize
                else:
                    p[0] = round(p[0]/gridSize)*gridSize
        verts = np.asarray(verts)
        primVor = Voronoi(verts)
        self.streetNodes = [StreetNode(i) for i in primVor.vertices]
        self.blocks = [Block(self) for i in primVor.regions]
        for i in range(len(self.blocks)):
            b = self.blocks[i]
            b.col = self.landColor
            for j in range(len(primVor.regions[i])):
                index = primVor.regions[i][j]
                if index >= 0:
                    b.verts.append(self.streetNodes[index])
        self.streets = []
        for i in range(len(primVor.ridge_vertices)):
            in0 = primVor.ridge_vertices[i][0]
            in1 = primVor.ridge_vertices[i][1]
            newStreet = Street(self.streetNodes[in0],self.streetNodes[in1],3)
            self.streets.append(newStreet)
        for i in self.blocks:
            for j in self.blocks:
                if i.sharedNeighbors(j) > 0:
                    i.neighborize(j)
        for k in self.blocks:
            k.centroid()
        self.neighborPts = []
        for n in self.node.neighbors:
            xReal = n.x - self.node.x
            yReal = n.y - self.node.y
            realDist = n.dist(self.node)
            optimalDist = math.sqrt((self.xDim**2)+(self.yDim**2))
            distMultiplier = optimalDist/realDist
            dx = xReal*distMultiplier
            dy = yReal*distMultiplier
            buffer = 16
            xx = clamp((self.xDim/2) + dx,buffer,self.xDim-buffer)
            yy = clamp((self.yDim/2) + dy,buffer,self.yDim-buffer)
            s = StreetNode([xx,yy])
            if n.river != None:
                s.type = "river"
                #s.drawCol = n.myMap.biomeColors[n.biome]
                s.drawCol = n.biomeColor
            if n.bodyWater != None:
                s.type = "water"
                s.drawCol = self.waterColor
            elif self.node.elevation > n.elevation:
                s.drawCol = n.biomeColor
                #s.drawCol = n.myMap.biomeColors[n.biome]
                #s.drawCol = colAvg(n.myMap.biomeColors[n.biome],self.landColor)
            else:
                s.drawCol = self.landColor
            if len(n.roads) > 0:
                s.road = True
            s.node = n
            self.neighborPts.append(s)
        for b in self.blocks:
            d = 10000;
            nearestNeighbor = self.neighborPts[0]
            for p in self.neighborPts:
                if b.blockDist(p.x,p.y) < d:
                    d = b.blockDist(p.x,p.y)
                    nearestNeighbor = p
            if b.blockDist(self.xDim/2,self.yDim/2) > d*0.75:
                b.type = nearestNeighbor.type
                b.col = nearestNeighbor.drawCol
        if self.node.river != None:
            riverPrevious = self.node.riverPrevious()
            riverNext = self.node.riverNext()
            for p in self.neighborPts:
                if p.node == riverPrevious or p.node == riverNext or p.type == "water":
                    self.buildRiver(p,self.nearestBlock(self.xDim/2,self.yDim/2))
        for b in self.blocks:
            for f in self.blocks:
                if (b.sharedNeighbors(f) >= 2 and
                    b.type == "water"):
                    f.col = self.waterColor
        if self.node.watery() == 1:
            for f in self.blocks:
                f.col = self.waterColor
        else:
            if self.nearestBlock(self.x,self.y).col == self.waterColor:
                newCenter = self.nearestSolidBlock(self.x,self.y)
                newCenter.centroid()
                self.x = newCenter.x
                self.y = newCenter.y
        for f in self.blocks:
            if f.col == self.waterColor:
                f.type = "water"
        self.avgColors(6)
    def nearestBlock(self,xx,yy):
        d = 100000
        a = self.blocks[0]
        for k in self.blocks:
            dd = k.blockDist(xx,yy)
            if dd < d and len(k.neighbors) > 0:
                d = dd
                a = k
        return a
    def nearestSolidBlock(self,xx,yy):
        d = 100000
        a = self.blocks[0]
        for k in self.blocks:
            dd = k.blockDist(xx,yy)
            if dd < d and len(k.neighbors) > 0 and k.col != self.waterColor:
                d = dd
                a = k
        return a
    def buildRiver(self,n0,b1):
        current = self.nearestBlock(n0.x,n0.y)
        bx = b1.centroid()[0]
        by = b1.centroid()[1]
        d = 10000
        while current != b1:
            choice = current
            for n in current.neighbors:
                if n.blockDist(bx,by) < d and len(n.neighbors) > 1:
                    choice = n
                    d = n.blockDist(bx,by)
            if choice == current:
                current = b1
                current.type = "water"
                current.color = self.waterColor
            else:
                current = choice
                current.type = "water"
                current.color = self.waterColor
    def avgColors(self,count=1):
        for u in range(count):
            for i in self.blocks:
                if i.col not in [self.waterColor,self.streetColor,"black"]:
                    c0 = i.col[0]
                    c1 = i.col[1]
                    c2 = i.col[2]
                    pp = 1
                    for j in i.neighbors:
                        if j.col not in [self.waterColor,self.streetColor,"black"]:
                            c0 += j.col[0]
                            c1 += j.col[1]
                            c2 += j.col[2]
                            pp += 1
                    c0 = math.floor(c0/pp)
                    c1 = math.floor(c1/pp)
                    c2 = math.floor(c2/pp)
                    i.col = (c0,c1,c2)
    def addBuilding(self,drawer,roadCenter):
        chosenX = math.floor(random.uniform(self.xDim/3,self.xDim/1.5))
        chosenY = math.floor(random.uniform(self.yDim/3,self.yDim/1.5))
        chosenBlock = self.nearestSolidBlock(chosenX,chosenY)
        if len(roadCenter) > 0:
            dCol = Tools.streetColor
            drawer.line([roadCenter,(chosenBlock.x,chosenBlock.y)],dCol,math.floor(self.roadRadius*2))
        chosenBlock.drawBuilding(drawer)
    def addFarm(self,drawer,roadCenter):
        polygonSides = random.choice([3,4,4,4,5,5,6])
        polygonRadius = math.floor(random.uniform(self.xDim/6,self.xDim/4))
        polygonAngle = random.randint(0,359);
        #polygonLocation = random.choice([(self.xDim/3,self.yDim/2),(self.xDim/1.5,self.yDim/2),(self.xDim/2,self.yDim/3),(self.xDim/2,self.yDim/1.5)])
        polygonCorners = []
        for p in range(polygonSides):
            lengthDeviation = random.uniform(0.65,1.2)
            corner_x = self.x + lengthDirX(polygonRadius*lengthDeviation,polygonAngle+((360/polygonSides)*p))
            corner_y = self.y + lengthDirY(polygonRadius*lengthDeviation,polygonAngle+((360/polygonSides)*p))
            polygonCorners.append((math.floor(corner_x),math.floor(corner_y)))
        chosenCorner = polygonCorners[0]
        chosenBlock = self.nearestSolidBlock(chosenCorner[0],chosenCorner[1])
        if len(roadCenter) > 0:
            dCol = Tools.streetColor
            drawer.line([roadCenter,(chosenBlock.x,chosenBlock.y)],dCol,math.floor(self.roadRadius*2))
        dCol = Tools.streetColor
        bCol = Tools.farmColor
        drawer.polygon(polygonCorners,bCol,dCol)
        polygonCorners.append(polygonCorners[0])
        drawer.line(polygonCorners,dCol,3)
        chosenBlock.drawBuilding(drawer)
    def addFort(self,drawer):
        polygonSides = random.choice([3,4,4,4,4,5,5,5,6,6,7,8])
        polygonRadius = math.floor(random.uniform(self.xDim/15,self.xDim/11))
        polygonAngle = random.randint(0,359);
        #polygonLocation = random.choice([(self.xDim/3,self.yDim/2),(self.xDim/1.5,self.yDim/2),(self.xDim/2,self.yDim/3),(self.xDim/2,self.yDim/1.5)])
        polygonCorners = []
        for p in range(polygonSides):
            corner_x = self.x + lengthDirX(polygonRadius,polygonAngle+((360/polygonSides)*p))
            corner_y = self.y + lengthDirY(polygonRadius,polygonAngle+((360/polygonSides)*p))
            polygonCorners.append((math.floor(corner_x),math.floor(corner_y)))
        dCol = Tools.streetColor
        bCol = Tools.buildingColor
        drawer.polygon(polygonCorners,bCol,dCol)
        polygonCorners.append(polygonCorners[0])
        drawer.line(polygonCorners,dCol,3)
    def drawRoads(self,image):
        drawer = ImageDraw.Draw(image)
        dCol = Tools.streetColor
        numRoads = len([i for i in self.neighborPts if i.road == True])
        roadCenter = ()
        if len(self.node.roads) > 0:
            if self.node.city != None or numRoads > 0:
                roadCenter = (self.x,self.y)
                for neighbor in self.neighborPts:
                    if neighbor.road == True:
                        sharedRoad = False
                        if neighbor.node in self.node.roads:
                            sharedRoad = True
                        if sharedRoad == True:
                            drawer.line([(neighbor.x,neighbor.y),(self.x,self.y)],dCol,math.floor(self.roadRadius*2))
            """
            else:
                for neighbor in self.neighborPts:
                    if neighbor.road == True:
                        for otherNeighbor in self.neighborPts:
                            if otherNeighbor.x != neighbor.x and otherNeighbor.y != neighbor.y and otherNeighbor.road == True:
                                drawer.line([(neighbor.x,neighbor.y),(otherNeighbor.x,otherNeighbor.y)],dCol,math.floor(self.roadRadius*2))
                                roadCenter = ((neighbor.x+otherNeighbor.x)/2,(neighbor.y+otherNeighbor.y)/2)
            """
        for k in self.blocks:
            if k.blockDist(self.x,self.y) < self.radius and k.col != Tools.waterColor:
                k.drawSelf(drawer)
        for k in self.blocks:
            if k.blockDist(self.x,self.y) < self.radius and k.col != Tools.waterColor:
                k.drawRoads(drawer)
        nodeStructure = self.node.structure()
        if self.node.city == None and nodeStructure != None:
            if nodeStructure in ["farm","mill"]:
                self.addFarm(drawer,roadCenter)
            if nodeStructure == "fort":
                self.addFort(drawer)
            if nodeStructure in ["inn","brothel","workshop","mine","fishery","port"]:
                self.addBuilding(drawer,roadCenter)
        scl = 1/BORDERSCALE
        h = self.yDim/scl
        w = self.xDim/scl
        dCol = (0,0,0)
        drawer.rectangle([(0,0),(self.xDim,h)],dCol,dCol)
        drawer.rectangle([(0,0),(w,self.yDim)],dCol,dCol)
        drawer.rectangle([(self.xDim-w,0),(self.xDim,self.yDim)],dCol,dCol)
        drawer.rectangle([(0,self.yDim-h),(self.xDim,self.yDim)],dCol,dCol)
    def drawSelf(self,drawer):
        drawer.rectangle([(0,0),(self.xDim,self.yDim)],self.landColor,self.landColor)
        for k in self.blocks:
            k.drawSelf(drawer)
            