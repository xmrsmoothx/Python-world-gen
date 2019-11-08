# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 18:02:07 2019

@author: Bri
"""

import random
import math

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

class Face:
    def __init__(self,c,n=None,p=None,x=128):
        self.culture = c
        self.node = n
        self.pop = p
        self.xDim = x
        self.bgc = (255,0,255)
    def generateCultureFace(self):
        random.seed(self.culture.name)
        inf = self.culture.value.influences.envInfluences
        hotness = inf["temperature"] # Typically between 0.200 and 0.800
        latitude = inf["latitude"] # Typically between 0.100 and 1.100
        humidity = inf["rainfall"] # Typically between 0.000 and 0.500
        # Hot climates produce thick lips and more skeletal builds.
        # Cold climates produce rounder, wider builds and long straight hair.
        # Sunny climates produce dark skin, hair, and eyes, and short curly hair.
        # Dark climates produce light skin, hair, and eyes.
        # Dry climates produce thin noses.
        # Humid climates produce wide noses and thick lips.
        self.headheight = clamp(0.9*self.xDim*random.uniform(0.7,1.1)*(math.sqrt(clamp(hotness,0.4,0.8))),self.xDim*0.6,self.xDim)
        self.headwidth = clamp(0.9*self.xDim*random.uniform(0.7,1.1)*(1-(hotness**2)),self.xDim/2.5,self.xDim*0.9)
        self.neckwidth = self.headwidth*0.6*random.uniform(0.8,1.2)
        self.eyesize = clamp(0.2*self.xDim*random.uniform(0.8,1.2)*(1-(0.5*latitude)),self.xDim/16,self.xDim/10)
        self.nosewidth = 0.15*self.xDim*random.uniform(0.7,1)*(math.sqrt(clamp(humidity,0.1,0.5)))
        self.lipthickness = clamp(0.2*self.xDim*random.uniform(0.8,1.2)*clamp(humidity,0.1,0.5)*clamp(hotness,0.1,0.7),self.xDim/24,self.xDim/7)
        self.lipwidth = self.lipthickness*3*random.uniform(0.7,1.3)
        self.melanin = clamp(254-(210*random.uniform(0.9,1.1)*(1.1-latitude)),42,230)
        self.skincolor = (math.floor(random.uniform(0.9,1.1)*14),clamp(math.floor(random.uniform(0.8,1.3)*(20+self.melanin)),24,200),clamp(math.floor(random.uniform(1,1.1)*(255-self.melanin)),32,250))
        self.shadowcolor = (self.skincolor[0],self.skincolor[1],math.floor(self.skincolor[2]*0.65))
        self.haircolor = (math.floor(random.uniform(0.8,1.4)*24),math.floor(random.uniform(0.9,1.1)*50),math.floor(random.uniform(0.9,1.1)*self.melanin))
        self.eyecolor = (math.floor(random.uniform(0.1,1)*240),math.floor(random.uniform(0.8,1.4)*180),250-(math.floor(random.uniform(0.6,1)*self.melanin)))
        self.hairlength = 0.7*self.xDim*random.uniform(0.8,1)*(1-hotness)
    def drawSelf(self,drawer):
        drawer.rectangle([(0,0),(self.xDim,self.xDim)],self.bgc,self.bgc)
        hmargin = (self.xDim-self.headwidth)/2
        vmargin = (self.xDim-self.headheight)/2
        drawer.rectangle([((self.xDim/2)-(self.neckwidth/2),vmargin+self.headheight/1.5),((self.xDim/2)+(self.neckwidth/2),self.xDim)],self.skincolor,self.skincolor)
        drawer.polygon([(0,self.xDim),(self.xDim,self.xDim),(self.xDim/2,self.xDim-vmargin)],self.skincolor,self.skincolor)
        drawer.ellipse([(hmargin,vmargin),(self.xDim-hmargin,vmargin+self.headheight/1.5)],self.skincolor)
        drawer.rectangle([(hmargin,vmargin+self.headheight/3),(self.xDim-hmargin,vmargin+self.headheight/1.5)],self.skincolor,self.skincolor)
        drawer.polygon([(hmargin,vmargin+self.headheight/1.5),(self.xDim-hmargin,vmargin+self.headheight/1.5),(self.xDim/2,self.xDim-vmargin)],self.skincolor,self.skincolor)
        drawer.line([(hmargin,vmargin+self.headheight/1.5),(self.xDim/2,self.xDim-vmargin),(self.xDim-hmargin,vmargin+self.headheight/1.5)],self.shadowcolor,1)
        scleracolor = (40,20,250)
        pupilcolor = (0,0,0)
        drawer.ellipse([((self.xDim/2)-(self.headwidth/3.5)-(self.eyesize/2),(self.xDim/2)-(self.eyesize/4.5)),
                        ((self.xDim/2)-(self.headwidth/3.5)+(self.eyesize/2),(self.xDim/2)+(self.eyesize/4.5))],scleracolor)
        drawer.ellipse([((self.xDim/2)+(self.headwidth/3.5)-(self.eyesize/2),(self.xDim/2)-(self.eyesize/4.5)),
                        ((self.xDim/2)+(self.headwidth/3.5)+(self.eyesize/2),(self.xDim/2)+(self.eyesize/4.5))],scleracolor)
        drawCircle(drawer,(self.xDim/2)-(self.headwidth/3.5),(self.xDim/2),self.eyesize/5,self.eyecolor)
        drawCircle(drawer,(self.xDim/2)+(self.headwidth/3.5),(self.xDim/2),self.eyesize/5,self.eyecolor)
        drawCircle(drawer,(self.xDim/2)-(self.headwidth/3.5),(self.xDim/2),self.eyesize/8,pupilcolor)
        drawCircle(drawer,(self.xDim/2)+(self.headwidth/3.5),(self.xDim/2),self.eyesize/8,pupilcolor)
        drawer.polygon([(self.xDim/2,self.xDim/2.1),(self.xDim/2,self.xDim/1.7),((self.xDim/2)+(self.nosewidth),self.xDim/1.8)],self.shadowcolor,self.shadowcolor)
        drawer.line([(self.xDim/2,self.xDim/1.7),((self.xDim/2)-(self.nosewidth),self.xDim/1.8)],self.shadowcolor,1)
        drawer.ellipse([((self.xDim/2)-(self.lipwidth/2)),((self.xDim/1.5)-(self.lipthickness/2)),
                        ((self.xDim/2)+(self.lipwidth/2)),((self.xDim/1.5)+(self.lipthickness/2))],self.shadowcolor)
        drawer.line([((self.xDim/2)-(self.lipwidth/2),self.xDim/1.5),((self.xDim/2)+(self.lipwidth/2),self.xDim/1.5)],self.skincolor,1)
        
        
        
        
        
        
        