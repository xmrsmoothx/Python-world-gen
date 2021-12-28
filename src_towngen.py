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
BORDERSCALE = 1/10
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
    def dist(self,x,y):
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
        return [self.x,self.y]
    def ccw(self,p0,p1,c):
        # Return 1 if p1 is counterclockwise from p0 with c as center; otherwise, return 0
        cx = c[0]
        cy = c[1]
        v0 = (p0.x-cx,p0.y-cy)
        v1 = (p1.x-cx,p1.y-cy)
        d0 = enorm(v0[0],v0[1])
        d1 = enorm(v1[0],v1[1])
        dp = (v0[0]*v1[0]) + (v0[1]*v1[1])
        if d0*d1 == 0:
            return 0
        q = clamp(dp/exclus(d0*d1),-1,1)
        ang = math.acos(q)
        if ang >= 0:
            return 1
        else:
            return 0
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
        n = len(self.verts)
        for i in range(n):
            for j in range(n):
                aa = j
                bb = (j+1) % n
                if self.ccw(self.verts[aa],self.verts[bb],self.centroid()) == 0:
                    # bb not ccw from aa
                    t = self.verts[aa]
                    self.verts[aa] = self.verts[bb]
                    self.verts[bb] = t
        self.verts = list(reversed(self.verts))
    def polyArea(self):
        if len(self.verts) < 1:
            self.area = 0
            return 0
        self.orderccw()
        area = 0
        q = self.verts[-1]
        for p in self.verts:
            area += (p.x * q.y) - (p.y * q.x)
            q = p
        area = abs(area/2)
        self.area = area
        return area
    def drawSelf(self,drawer):
        dCol = self.col
        if len(self.verts) > 1:
            vts = [(p.x,p.y) for p in self.verts]
            drawer.polygon(vts,dCol,dCol)

class Town:
    def __init__(self,n,m,nom):
        self.myMap = m
        self.node = n
        self.name = nom
        self.mapName = "./generated/town_"+self.name+".gif"
        self.xDim = XDIM
        self.yDim = XDIM
        self.landColor = self.myMap.biomeColors[self.node.biome]
        self.streetColor = Tools.streetColor
        self.waterColor = Tools.waterColor
        cnt = 512
        sd = (self.node.x*73)+(self.node.y*37)
        random.seed(sd)
        buffer = 32
        verts = [[random.randint(buffer,self.xDim-buffer),random.randint(buffer,self.yDim-buffer)] for i in range(cnt)]
        verts.append([self.xDim/2,self.yDim/2])
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
        if self.node.city != None:
            self.population = self.node.city.population
            self.wateryNeighbors = [u for u in self.node.neighbors if u.watery() == 1]
            self.radius = math.floor(28+(math.sqrt(self.population)))*(1+(len(self.wateryNeighbors)/3))
            if self.node.river != None:
                self.radius = self.radius+36
            polygonSides = 5 + (sd % 7)
            centerX = XDIM/2;
            centerY = XDIM/2;
            self.polygon = [[centerX+(math.sin(6.283*u/polygonSides)*self.radius)+((u*sd) % 31),centerY+(math.cos(6.283*u/polygonSides)*self.radius)+((u*sd) % 29)] for u in range(math.floor(polygonSides))]
        else:
            self.population = 0
            self.radius = 0
        self.x = self.xDim/2
        self.y = self.yDim/2
        self.neighborPts = []
        for n in self.node.neighbors:
            d = n.dist(self.node)
            scl = self.xDim/d
            dx = (n.x-self.node.x)*scl
            dy = (n.y-self.node.y)*scl
            while (abs(dx) < self.xDim/2 or abs(dy) < self.yDim/2) and n.river != None and n.bodyWater == None:
                dx *= 1.1
                dy *= 1.1
            buffer = 16
            xx = clamp(self.x + dx,buffer,self.xDim-buffer)
            yy = clamp(self.y + dy,buffer,self.yDim-buffer)
            s = StreetNode([xx,yy])
            if n.river != None:
                s.type = "river"
                s.drawCol = self.waterColor
            if n.bodyWater != None:
                s.type = "water"
                s.drawCol = self.waterColor
            else:
                s.drawCol = n.myMap.biomeColors[n.biome]
                #s.drawCol = colAvg(n.myMap.biomeColors[n.biome],self.landColor)
            if len(n.roads) > 0:
                s.road = True
            self.neighborPts.append(s)
        for p in self.neighborPts:
            if p.type != "river":
                for k in range(len(self.blocks)):
                    kk = self.blocks[k].centroid()
                    if (self.blocks[k].blockDist(p.x,p.y) < self.x*0.9):
                        self.blocks[k].type = p.type
                        self.blocks[k].col = p.drawCol
            if p.type == "river"  or p.type == "water":
                if self.node.river != None:
                    self.buildRiver(p,self.nearestBlock(self.x,self.y))
        for b in self.blocks:
            for f in self.blocks:
                if (b.sharedNeighbors(f) >= 2 and
                    b.type == "water"):
                    f.col = self.waterColor
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
    def buildRiver(self,n0,b1):
        #probably scrap this and remake
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
    def drawRoads(self,image):
        drawer = ImageDraw.Draw(image)
        dCol = Tools.streetColor
        bCol = Tools.buildingColor
        centerX = XDIM/2
        centerY = XDIM/2
        if self.node.city != None:
            matPath = mpltPath.Path(self.polygon)
            for xx in range(math.floor(XDIM-(XDIM*BORDERSCALE*2))):
                for yy in range(math.floor(XDIM-(XDIM*BORDERSCALE*2))):
                    x = xx+(XDIM*BORDERSCALE)
                    y = yy+(XDIM*BORDERSCALE)
                    if matPath.contains_point((x,y)):
                        if (image.getpixel((x,y))) not in [(0,0,0),self.waterColor]:
                            if xStreets(x,self.name) or yStreets(y,self.name):
                                drawer.point((x,y),dCol)
                            else:
                                drawer.point((x,y),bCol)
        roadRadius = 1.5
        for neighbor in self.neighborPts:
            if neighbor.road == True:
                drawCircle(drawer,centerX,centerY,roadRadius,dCol)
                drawer.line([(neighbor.x,neighbor.y),(centerX,centerY)],dCol,math.floor(roadRadius*2))
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
            