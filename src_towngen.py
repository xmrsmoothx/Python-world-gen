# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 14:52:37 2017

@author: Bri
"""

import math
import random

def clamp(x,minimum,maximum):
    if x < minimum:
        return minimum
    elif x > maximum:
        return maximum
    else:
        return x

def drawCircle(drawer,x,y,radius,color):
    x1 = x-radius
    x2 = x+radius
    y1 = y-radius
    y2 = y+radius
    drawer.ellipse([(x1,y1),(x2,y2)],color)
    
def lengthDirX(length, angle):
  radian_angle = math.radians(angle)
  return length * math.cos(radian_angle)

def lengthDirY(length, angle):
  radian_angle = math.radians(angle)
  return length * math.sin(radian_angle)

def A(dx, dy):
  return math.degrees( math.atan2(dy, dx) )

def drawTrapezoid(drawer,x1,y1,x2,y2,r1,r2,color):
    directAngle = A(x2-x1,y2-y1)
    pAngle = directAngle-90
    pAngle2 = pAngle-180
    p1 = (x1+lengthDirX(r1,pAngle),y1+lengthDirY(r1,pAngle))
    p2 = (x1+lengthDirX(r1,pAngle2),y1+lengthDirY(r1,pAngle2))
    p3 = (x2+lengthDirX(r2,pAngle2),y2+lengthDirY(r2,pAngle2))
    p4 = (x2+lengthDirX(r2,pAngle),y2+lengthDirY(r2,pAngle))
    drawer.polygon([p1,p2,p3,p4],color,color)

def getKey(n):
    return n.key

def sd(node):
    seed = 0
    seed += node.x*0.1
    seed += node.y*0.25
    seed += node.elevation*8
    seed += (node.x*node.y) % 17
    return seed

def ptDist(n1,n2,dt):
    d = math.sqrt(((n2.x-n1.x)**2)+((n2.y-n1.y)**2))
    t = dt/d
    xx = ((1-t)*n1.x) + (t*n2.x)
    yy = ((1-t)*n1.y) + (t*n2.y)
    n = StreetNode(xx,yy,n1.town)
    return n

def binRoll(n,d):
    if math.floor(n.x) % 2 == 0:
        s = (n.x*89*d) + (n.y*73*d)
    else:
        s = (n.x*37*d) + (n.y*53*d)
    s = (s % 31) + d
    if math.floor(s) % math.floor(d) == 0:
        r = 1
    else:
        r = 0
    return r

def seedAmount(n,d,minimum,maximum):
    d = d + 1
    if math.floor(n.x) % 2 == 0:
        s = (n.x*149*d) + (n.y*17*d)
    else:
        s = (n.x*73*d) + (n.y*29*d)
    ct = maximum-minimum
    val = (s % ct) + minimum
    return val

class StreetNode:
    def __init__(self,xx,yy,t,k=0):
        self.town = t
        self.x = xx
        self.y = yy
        self.key = k
        self.neighbors = []
        self.outRoad = 0
        self.outRiver = 0
    def link(self,n):
        if n not in self.neighbors:
            self.neighbors.append(n)
            n.neighbors.append(n)
            e = StreetEdge(self,n,self.town)
    def isLinked(self,n):
        if n in self.neighbors:
            return 1
        else:
            return 0
    def coords(self):
        return (self.x,self.y)
    def perturb(self):
        pDist = math.ceil(self.town.xDim/64)
        xDrift = (self.x*71+self.y*83) % pDist
        if math.floor(xDrift) % 2 == 0:
            xDrift *= -1
        yDrift = (self.x*83+self.y*71) % pDist
        if math.floor(yDrift) % 2 == 0:
            yDrift *= -1
        self.x = self.x+xDrift
        self.y = self.y+yDrift
    def drawSelf(self,drawer):
        drawCircle(drawer,self.x,self.y,3,(0,0,0))

class StreetEdge:
    def __init__(self,n1,n2,t):
        self.node1 = n1
        self.node2 = n2
        self.town = t
        self.town.edges.append(self)
    def drawSelf(self,drawer):
        dCol = (0,0,0)
        r1 = 1
        r2 = 1
        drawTrapezoid(drawer,self.node1.x,self.node1.y,self.node2.x,self.node2.y,r1,r2,dCol)
        
class TownHex:
    def __init__(self,size,t,hexid=0):
        self.town = t
        self.verts = self.town.verts
        self.pts = []
        self.innerPts = []
        self.town.d = self.town.d + seedAmount(self.town.center,size+3,1,3)
        for i in range(self.verts):
            dist = (size+1+self.town.d)*((self.town.xDim*0.8)/self.town.maxHexSize)
            n = ptDist(self.town.center,self.town.nbrs[i],dist)
            self.pts.append(n)
            r = binRoll(n,3)
            if size > 0:
                parent = self.town.hexes[size-1]
                self.innerPts.append(parent.pts[i])
                if r == 0:
                    self.pts[i].link(self.innerPts[i])
            r = binRoll(n,3.3)
            if i != 0:
                if r == 0:
                    self.pts[i].link(self.pts[i-1])
            if i == self.verts-1:
                if r == 0:
                    self.pts[i].link(self.pts[0])
        for q in self.pts:
            q.x = (self.town.center.x+q.x)/2
            q.y = (self.town.center.y+q.y)/2
        for q in self.pts:
            q.perturb()
                
    
class Town:
    def __init__(self,n,m,nom=""):
        self.node = n
        self.myMap = m
        self.mapName = "town_"+nom+".gif"
        self.xDim = 512
        self.yDim = 512
        xx = self.xDim/2
        yy = self.yDim/2
        self.center = StreetNode(xx,yy,self)
        self.nbrs = []
        self.d = 0
        for k in self.node.neighbors:
            scale = 16
            dx = (k.x-self.node.x)*scale
            dy = (k.y-self.node.y)*scale
            x1 = self.center.x+dx
            y1 = self.center.y+dy
            key = math.atan2(y1-yy,x1-xx)
            nd = StreetNode(x1,y1,self,key)
            if k in self.node.roads:
                nd.outRoad = 1
            if k.river == self.node.river:
                nd.outRiver = 1
            self.nbrs.append(nd)
        self.nbrs = sorted(self.nbrs,key=lambda streetnode: streetnode.key)
        self.verts = len(self.nbrs)
        self.seed = sd(self.node)
        self.population = self.node.city.population
        self.maxHexSize = 16
        p = 0
        self.hexes = []
        self.edges = []
        while p < clamp(math.log2(self.population)-4,2,self.maxHexSize):
            newHex = TownHex(p,self,hexid=p)
            self.hexes.append(newHex)
            p+=1
    def drawSelf(self,drawer):
        for e in self.edges:
            e.drawSelf(drawer)