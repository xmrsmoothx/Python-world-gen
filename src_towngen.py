# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 14:52:37 2017

@author: Bri
"""

import numpy as np
import random
import math
from scipy.spatial import Voronoi

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
        return d
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
        v0x = p0.x-cx
        v0y = p0.y-cy
        v1x = p1.x-cx
        v1y = p1.y-cy
        if v0y*v1x > v0x*v1y:
            return 1
        else:
            return 0
    def blockDist(self,x,y):
        self.centroid()
        dx = abs(self.x-x)
        dy = abs(self.y-y)
        dist = math.sqrt((dx**2) + (dy**2))
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
                    t = self.verts[aa]
                    self.verts[aa] = self.verts[bb]
                    self.verts[bb] = t
    def polyArea(self):
        self.orderccw()
        n = len(self.verts)
        area = 0
        for i in range(n):
            j = (i+1) % n
            aa = self.verts[i]
            bb = self.verts[j]
            yy = (aa.y+bb.y)/2
            xx = (aa.x-bb.x)
            zz = xx*yy
            area += zz
        area = abs(area)
        self.area = area
        return area
    def subdivide(self,size=1028,level=0):
        if level > 5:
            return -1
        a = self.polyArea()
        if self.area > size and len(self.verts) > 3:
            random.seed(self.area)
            sub = 0
            div0 = None
            div1 = None
            pts = []
            pts1 = []
            while sub == 0:
                base0 = random.choice(self.verts)
                nbr0 = random.choice(base0.neighbors)
                while nbr0 not in self.verts:
                    nbr0 = random.choice(base0.neighbors)
                div0 = StreetNode(base0.midpt(nbr0))
                s = Street(div0,nbr0)
                s = Street(div0,base0)
                base1 = random.choice(self.verts)
                nbr1 = random.choice(base1.neighbors)
                while nbr1 not in self.verts:
                        nbr1 = random.choice(base1.neighbors)
                div1 = StreetNode(base1.midpt(nbr1))
                while div1 == div0:
                    base1 = random.choice(self.verts)
                    nbr1 = random.choice(base1.neighbors)
                    while nbr1 not in self.verts:
                        nbr1 = random.choice(base1.neighbors)
                    div1 = StreetNode(base1.midpt(nbr1))
                s = Street(div1,nbr1)
                s = Street(div1,base1)
                c = div0.midpt(div1)
                pts = []
                for p in self.verts:
                    if self.ccw(div1,p,c) == 0 and self.ccw(div0,p,c) == 1:
                        pts.append(p)
                pts1 = []
                for p in self.verts:
                    if p not in pts:
                        pts1.append(p)
                l0 = len(pts)+2
                l1 = len(pts1)+2
                if l0 > 1 and l1 > 1:
                    sub = 1
            ww = random.randint(2,2)
            q = Street(div0,div1,ww)
            s = Block(self.myTown)
            s.verts = [div1,div0]
            s.verts.extend(pts)
            s.subdivide(size,level+1)
            self.subblocks.append(s)
            s1 = Block(self.myTown)
            s1.verts = [div1,div0]
            s1.verts.extend(pts1)
            s1.subdivide(size,level+1)
            self.subblocks.append(s1)
            self.substreets.append(q)
    def drawSelf(self,drawer):
        dCol = self.col
        if len(self.verts) > 1:
            vts = [(p.x,p.y) for p in self.verts]
            drawer.polygon(vts,dCol,dCol)
        for k in self.subblocks:
            k.drawSelf(drawer)
        for k in self.substreets:
            k.drawSelf(drawer,self.myTown.streetColor)

class Town:
    def __init__(self,n,m,nom):
        self.myMap = m
        self.node = n
        self.name = nom
        self.mapName = "town_"+self.name+".gif"
        self.xDim = 720
        self.yDim = 720
        self.landColor = self.myMap.biomeColors[self.node.biome]
        self.streetColor = (16,128,76)
        self.waterColor = (142,78,64)
        cnt = 256
        sd = (self.node.x*73)+(self.node.y*37)
        random.seed(sd)
        verts = [[random.randint(0,self.xDim),random.randint(0,self.yDim)] for i in range(cnt)]
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
        if self.node.city != None:
            self.population = self.node.city.population
            self.radius = 32+(math.sqrt(self.population)*2)
        else:
            self.population = 0
            self.radius = -8
        self.x = self.xDim/2
        self.y = self.yDim/2
        self.neighborPts = []
        for n in self.node.neighbors:
            d = n.dist(self.node)
            scl = self.xDim/d
            dx = (n.x-self.node.x)*scl
            dy = (n.y-self.node.y)*scl
            xx = clamp(self.x + dx,0,self.xDim)
            yy = clamp(self.y + dy,0,self.yDim)
            s = StreetNode([xx,yy])
            if n.river != None:
                s.type = "river"
                s.drawCol = self.waterColor
            if n.bodyWater != None:
                s.type = "water"
                self.radius *= 1.15
                s.drawCol = self.waterColor
            else:
                s.drawCol = colAvg(n.biomeColor,self.landColor)
            self.neighborPts.append(s)
        for p in self.neighborPts:
            if p.type == "river":
                if self.node.river != None:
                    self.buildRiver(p,self.nearestBlock(self.x,self.y))
            else:
                for k in range(len(self.blocks)):
                    kk = self.blocks[k].centroid()
                    if (self.blocks[k].blockDist(p.x,p.y) < self.x):
                        self.blocks[k].type = p.type
                        self.blocks[k].col = p.drawCol
        self.outskirts = []
        for k in range(len(self.blocks)):
            if (self.blocks[k].blockDist(self.x,self.y) > self.radius
                or self.blocks[k].type == "water"):
                self.outskirts.append(self.blocks[k])
        self.cullStreets(0)
        for i in self.blocks:
            for j in self.blocks:
                if i.sharedNeighbors(j) > 0:
                    i.neighborize(j)
        for i in self.blocks:
            if i not in self.outskirts:
                i.col = "black"
                i.subdivide(800)
        self.avgColors(3)
    def nearestBlock(self,xx,yy):
        d = self.xDim
        a = self.blocks[0]
        for k in self.blocks:
            dd = k.blockDist(xx,yy)
            if dd < d:
                d = dd
                a = k
        return a
    def buildRiver(self,n0,b1):
        curr = self.nearestBlock(n0.x,n0.y)
        bx = b1.centroid()[0]
        by = b1.centroid()[1]
        while curr != b1:
            curr.type = "water"
            curr.col = (142,78,64)
            t = curr
            d = self.xDim
            for f in self.blocks:
                if (curr.sharedNeighbors(f) >= 2 and f != curr
                and f.blockDist(bx,by) < d):
                    d = f.blockDist(bx,by)
                    t = f
            if t == curr:
                return -1
            else:
                curr = t
    def cullStreets(self,l):
        for k in self.streets:
            k.exists = 0
            q = 0
            f = 0
            for b in self.blocks:
                if b not in self.outskirts:
                    if k.nodes[0] in b.verts:
                        q = 1
                    if k.nodes[1] in b.verts:
                        f = 1
            if q == 1 and f == 1:
                k.exists = 1
            if k.length() < l:
                k.cull()
    def avgColors(self,count=1):
        for u in range(count):
            for i in self.blocks:
                if i.col not in [self.waterColor,"black"]:
                    c0 = i.col[0]
                    c1 = i.col[1]
                    c2 = i.col[2]
                    pp = 1
                    for j in i.neighbors:
                        if j.col not in [self.waterColor,"black"]:
                            c0 += j.col[0]
                            c1 += j.col[1]
                            c2 += j.col[2]
                            pp += 1
                    c0 = math.floor(c0/pp)
                    c1 = math.floor(c1/pp)
                    c2 = math.floor(c2/pp)
                    i.col = (c0,c1,c2)
    def drawSelf(self,drawer):
        drawer.rectangle([(0,0),(self.xDim,self.yDim)],self.landColor,self.landColor)
        for k in self.blocks:
            k.drawSelf(drawer)
        for s in self.streets:
            if s.exists == 1:
                s.drawSelf(drawer,self.streetColor)
        scl = 8
        h = self.yDim/scl
        w = self.xDim/scl
        dCol = (0,0,0)
        drawer.rectangle([(0,0),(self.xDim,h)],dCol,dCol)
        drawer.rectangle([(0,0),(w,self.yDim)],dCol,dCol)
        drawer.rectangle([(self.xDim-w,0),(self.xDim,self.yDim)],dCol,dCol)
        drawer.rectangle([(0,self.yDim-h),(self.xDim,self.yDim)],dCol,dCol)