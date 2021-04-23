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

def avg2(v1,v2):
    return ((v1+v1)/2)

class Face:
    def __init__(self,c,n=None,p=None,x=128):
        self.culture = c
        self.node = n
        self.pop = p
        self.xDim = x
        self.bgc = (255,0,255)
    def shiftVals(self):
        self.headheight = self.headheight*random.uniform(0.9,1.1)
        self.headwidth = self.headwidth*random.uniform(0.9,1.1)
        self.neckwidth = self.headwidth*0.6*random.uniform(0.9,1.1)
        self.eyesize = self.eyesize*random.uniform(0.9,1.1)
        self.nosewidth = self.nosewidth*random.uniform(0.9,1.1)
        self.lipthickness = self.lipthickness*random.uniform(0.9,1.1)
        self.lipwidth = self.lipthickness*3.7*random.uniform(0.9,1.1)
        self.melanin = math.floor(self.melanin*random.uniform(0.9,1.1))
        self.skincolor = (math.floor(self.skincolor[0]*random.uniform(0.9,1.1)),
                          math.floor(self.skincolor[1]*random.uniform(0.9,1.1)),
                          math.floor(self.skincolor[2]*random.uniform(0.9,1.1)))
        self.shadowcolor = (self.skincolor[0],self.skincolor[1],math.floor(self.skincolor[2]*0.65))
        self.haircolor = (math.floor(self.haircolor[0]*random.uniform(0.8,1.2)),
                          math.floor(self.haircolor[1]*random.uniform(0.8,1.2)),
                          math.floor(self.haircolor[2]*random.uniform(0.8,1.2)))
        self.eyecolor = (math.floor(self.eyecolor[0]*random.uniform(0.6,1.6)),
                          math.floor(self.eyecolor[1]*random.uniform(0.8,1.2)),
                          math.floor(self.eyecolor[2]*random.uniform(0.8,1.2)))
        self.hairlength = self.hairlength*random.uniform(0.85,1.15)
        self.hairthickness = self.hairthickness*random.uniform(0.8,1.2)
        self.partleft = ((self.hairlength/self.xDim)*2*60)*random.uniform(0.7,1.3)
        self.partright = ((self.hairlength/self.xDim)*2*60)*random.uniform(0.7,1.3)
        self.earsize = self.earsize*random.uniform(0.9,1.1)
    def clampVals(self):
        self.headheight = clamp(self.headheight,self.xDim*0.7,self.xDim*0.75)
        self.headwidth = clamp(self.headwidth,self.xDim*0.5,self.xDim*0.65)
        self.eyesize = clamp(self.eyesize,self.xDim/13,self.xDim/9)
        self.nosewidth = clamp(self.nosewidth,self.xDim*0.08,self.xDim*0.23)
        self.lipthickness = clamp(self.lipthickness,self.xDim/24,self.xDim/11)
        self.melanin = clamp(self.melanin,10,230)
        self.skincolor = (clamp(self.skincolor[0],11,17),clamp(self.skincolor[1],24,200),clamp(self.skincolor[2],43,247))
        self.haircolor = (clamp(self.haircolor[0],4,36),clamp(self.haircolor[1],160,240),clamp(self.haircolor[2],0,254))
        self.eyecolor = (clamp(self.eyecolor[0],5,240),clamp(self.eyecolor[1],144,252),clamp(self.eyecolor[2],5,250))
        self.hairlength = clamp(self.hairlength,0,self.xDim*0.7)
        self.hairthickness = clamp(self.hairthickness,1,7)
        self.earsize = clamp(self.earsize,6,10)
    def generateCultureFace(self):
        random.seed(self.culture.name)
        inf = self.culture.value.influences.envInfluences
        hotness = inf["temperature"] # Typically between 0.200 and 0.800
        latitude = inf["latitude"] # Typically between 0.100 and 1.100
        humidity = inf["rainfall"] # Typically between 0.000 and 0.500
        # Hot climates produce thick lips and more skeletal builds.
        # Cold climates produce rounder, wider builds and long thin straight hair.
        # Sunny climates produce dark skin, hair, and eyes, and short thick curly hair.
        # Dark climates produce light skin, hair, and eyes.
        # Dry climates produce thin noses.
        # Humid climates produce wide noses and thick lips.
        self.headheight = clamp(0.9*self.xDim*random.uniform(0.7,1.1)*(math.sqrt(clamp(hotness,0.4,0.8))),self.xDim*0.7,self.xDim*0.75)
        self.headwidth = clamp(0.9*self.xDim*random.uniform(0.7,1.1)*(1-(hotness**2)),self.xDim*0.5,self.xDim*0.65)
        self.neckwidth = self.headwidth*0.6*random.uniform(0.8,1.2)
        self.eyesize = clamp(0.25*self.xDim*random.uniform(0.8,1.2)*(1-(0.5*latitude)),self.xDim/12,self.xDim/8)
        self.nosewidth = 0.12*self.xDim*random.uniform(0.7,1)*(math.sqrt(clamp(humidity,0.1,0.5)))
        self.lipthickness = clamp(0.23*self.xDim*random.uniform(0.8,1.2)*clamp(humidity,0.1,0.5)*clamp(hotness,0.1,0.7),self.xDim/20,self.xDim/10)
        self.lipwidth = self.lipthickness*3.7*random.uniform(0.7,1.3)
        self.melanin = clamp(254-(230*random.uniform(0.9,1.1)*(1.1-latitude)),10,230)
        self.skincolor = (math.floor(random.uniform(0.85,1.15)*14),clamp(math.floor(random.uniform(0.8,1.3)*(20+self.melanin)),24,200),clamp(math.floor(random.uniform(1,1.1)*(255-self.melanin)),43,247))
        self.shadowcolor = (self.skincolor[0],self.skincolor[1],math.floor(self.skincolor[2]*0.65))
        self.haircolor = (math.floor(random.uniform(-16,16)+20),math.floor(random.uniform(0.8,1.2)*200),clamp(math.floor(random.uniform(0.9,1.1)*(200-self.melanin)),0,254))
        self.eyecolor = (math.floor(random.uniform(0.01,1)*250),math.floor(random.uniform(0.8,1.4)*180),250-(math.floor(random.uniform(0.6,1)*self.melanin)))
        self.hairlength = math.floor(0.7*self.xDim*random.uniform(0.8,1)*(1-hotness))
        self.hairthickness = math.floor(1.8+(hotness*5.3))
        self.partleft = ((self.hairlength/self.xDim)*2*60)*random.uniform(0.7,1.3)
        self.partright = ((self.hairlength/self.xDim)*2*60)*random.uniform(0.7,1.3)
        self.earsize = (6+(5.3*humidity))*random.uniform(0.8,1.2)
        self.clampVals()
    def generateFace1(self,parent):
        self.headheight = parent.headheight
        self.headwidth = parent.headwidth
        self.neckwidth = parent.neckwidth
        self.eyesize = parent.eyesize
        self.nosewidth = parent.nosewidth
        self.lipthickness = parent.lipthickness
        self.lipwidth = parent.lipwidth
        self.melanin = parent.melanin
        self.skincolor = parent.skincolor
        self.shadowcolor = parent.shadowcolor
        self.haircolor = parent.haircolor
        self.eyecolor = parent.eyecolor
        self.hairlength = parent.hairlength
        self.hairthickness = parent.hairthickness
        self.earsize = parent.earsize
        self.shiftVals()
        self.clampVals()
    def generateFace2(self,parent1,parent2):
        random.seed(self.pop.name)
        self.headheight = avg2(parent1.headheight,parent2.headheight)
        self.headwidth = avg2(parent1.headwidth,parent2.headwidth)
        self.neckwidth = avg2(parent1.neckwidth,parent2.neckwidth)
        self.eyesize = avg2(parent1.eyesize,parent2.eyesize)
        self.nosewidth = avg2(parent1.nosewidth,parent2.nosewidth)
        self.lipthickness = avg2(parent1.lipthickness,parent2.lipthickness)
        self.lipwidth = avg2(parent1.lipwidth,parent2.lipwidth)
        self.melanin = avg2(parent1.melanin,parent2.melanin)
        self.skincolor = (math.floor(avg2(parent1.skincolor[0],parent2.skincolor[0])),
                          math.floor(avg2(parent1.skincolor[1],parent2.skincolor[1])),
                          math.floor(avg2(parent1.skincolor[2],parent2.skincolor[2])))
        self.shadowcolor = (self.skincolor[0],self.skincolor[1],math.floor(self.skincolor[2]*0.65))
        self.haircolor = (math.floor(avg2(parent1.skincolor[0],parent2.skincolor[0])),
                          math.floor(avg2(parent1.skincolor[1],parent2.skincolor[1])),
                          math.floor(avg2(parent1.skincolor[2],parent2.skincolor[2])))
        self.eyecolor = (math.floor(avg2(parent1.eyecolor[0],parent2.eyecolor[0])),
                          math.floor(avg2(parent1.eyecolor[1],parent2.eyecolor[1])),
                          math.floor(avg2(parent1.eyecolor[2],parent2.eyecolor[2])))
        self.hairlength = avg2(parent1.hairlength,parent2.hairlength)
        self.hairthickness = avg2(parent1.hairthickness,parent2.hairthickness)
        self.earsize = avg2(parent1.earsize,parent2.earsize)
        self.shiftVals()
        self.clampVals()
    def drawSelf(self,drawer):
        drawer.rectangle([(0,0),(self.xDim,self.xDim)],self.bgc,self.bgc)
        hmargin = (self.xDim-self.headwidth)/2
        vmargin = (self.xDim-self.headheight)/2
        drawer.ellipse([(hmargin-self.hairthickness,vmargin-self.hairthickness),
                        (self.xDim-(hmargin-self.hairthickness),vmargin+self.headheight/1.5)],self.haircolor)
        drawer.rectangle([((self.xDim/2)-(self.neckwidth/2),vmargin+self.headheight/1.5),((self.xDim/2)+(self.neckwidth/2),self.xDim)],self.skincolor,self.skincolor)
        drawer.polygon([(0,self.xDim),(self.xDim,self.xDim),(self.xDim/2,self.xDim-vmargin)],self.skincolor,self.skincolor)
        drawer.ellipse([(hmargin,vmargin),(self.xDim-hmargin,vmargin+self.headheight/1.5)],self.skincolor)
        drawer.chord([(hmargin-1,vmargin-1),(self.xDim+1-hmargin,vmargin+1+self.headheight/1.5)],
                      240-self.partleft,300+self.partright,self.haircolor,self.haircolor)
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
        drawer.rectangle([((self.xDim/2)-(self.headwidth/2),(self.xDim/2)-(self.eyesize/1.5)),
                          ((self.xDim/2)+(self.headwidth/2),(self.xDim/2)-(self.eyesize/5))],self.skincolor,self.skincolor)
        drawer.ellipse([(hmargin+1-self.earsize,(self.xDim/2)+3-(self.earsize*2)),
                        (hmargin+3,(self.xDim/2)+self.earsize)],self.skincolor,self.skincolor)
        drawer.ellipse([(self.xDim-(hmargin+3),(self.xDim/2)+3-(self.earsize*2)),
                        ((self.xDim-(hmargin+1))+self.earsize,(self.xDim/2)+self.earsize)],self.skincolor,self.skincolor)
        
        
        
        
        
        
        