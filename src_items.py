# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 01:07:42 2018

@author: Bri
"""

import random
from PIL import Image, ImageFont, ImageDraw, ImageTk
from src_tools import *
from src_events import *
import string
import math

class Sigil:
    def __init__(self,drawer,c,col,a=random.choice([-90,90]),s=None,o=None):
        if s != None:
            random.seed(s)
        self.owner = o
        self.center = c
        self.components = []
        self.numComponents = 0
        self.seed = s
        maxComponents = 3
        minComponents = 0
        maxRadius = 70
        if self.owner != None:
            maxComponents = math.floor(self.owner.numComponents/2)
            maxRadius = math.floor(self.owner.radius/2.5)
        if maxRadius <= 3:
            return
        self.radius = random.randint(3,maxRadius)
        self.edgeCount = random.randint(-3,12)
        self.circle = False
        if self.edgeCount <= 2:
            self.edgeCount = random.randint(3,12)
            self.circle = True
        if maxComponents > 0:
            self.numComponents = random.randint(minComponents,maxComponents)
        for componentIndex in range(self.numComponents):
            maxCount = self.edgeCount
            componentCount = random.choice(factorsOf(maxCount))
            componentDist = random.randint(-math.floor(self.radius*0.4),math.floor(self.radius*0.33))
            self.components.append({"count":componentCount,"dist":componentDist})
        angleIncrement = 360/self.edgeCount
        currentAngle = a
        self.componentsOnEdges = False
        if random.random() > 0.5:
            self.componentsOnEdges = True
            currentAngle = currentAngle - (angleIncrement*0.5)
        self.star = False
        if random.random() > 0.45 and self.edgeCount >= 5:
            self.star = True
        self.lineThickness = random.randint(1,3)
        for vertexIndex in range(0,self.edgeCount):
            for component in self.components:
                componentIndex = (vertexIndex/self.edgeCount)*component["count"]
                if componentIndex-math.floor(componentIndex) < 0.002:
                    componentAngle = currentAngle
                    componentDist = math.floor(self.radius+component["dist"])
                    if self.componentsOnEdges == True:
                        componentAngle += angleIncrement*0.5
                    componentCenter = (self.center[0]+lengthDirX(componentDist,componentAngle),self.center[1]+lengthDirY(componentDist,componentAngle))
                    addSigil = Sigil(drawer,componentCenter,col,a=componentAngle,s=self.seed*componentDist,o=self)
            currentVertex = (self.center[0]+lengthDirX(self.radius,currentAngle),self.center[1]+lengthDirY(self.radius,currentAngle))
            nextVertexAngle = currentAngle+angleIncrement
            if self.star == True:
                nextVertexAngle = currentAngle+(angleIncrement*2)
            nextVertex = (self.center[0]+lengthDirX(self.radius,nextVertexAngle),self.center[1]+lengthDirY(self.radius,nextVertexAngle))
            if self.circle == False:
                drawer.line([currentVertex,nextVertex],col,self.lineThickness)
            currentAngle += angleIncrement
        if self.circle == True:
            drawCircle(drawer,self.center[0],self.center[1],self.radius,col,out=True)
                    

class Item:
    def __init__(self,k,c,nn="",f=None,s=None,i=3,cr=None):
        self.tt = "item"
        self.kind = k
        self.culture = c
        self.culture.items.append(self)
        self.gender = 0
        self.creationEvent = None
        self.destructionEvent = None
        self.condition = 1
        self.field = f
        self.subject = s
        self.importance = i
        self.location = None
        self.move(self.culture.origin)
        self.creator = cr
        self.owner = None
        if self.creator != None:
            self.move(self.creator.location)
            self.creator.inventory.append(self)
            self.owner = self.creator
        self.subkind = None
        self.decoration = None
        self.quality = random.uniform(0.1,0.3)
        self.quality = clamp(self.quality+(self.creator.skill*random.uniform(0.75,1.1)),0.05,1)
        self.importance = self.importance*(1+(self.quality/2))
        if self.kind == "piece":
            self.subkind = synonym("piece")
        subname = ""
        if self.subject != None:
            if self.subject.tt == "event":
                self.importance = clamp(self.importance+math.sqrt(self.subject.importance)/2,self.importance,self.importance*3)
                subname = self.subject.note()
            if self.subject.tt == "pop":
                self.importance = clamp(self.importance+math.sqrt(self.subject.importance)/2,self.importance,self.importance*3)
                subname = "The " + self.subject.nameFull()
            if self.subject.tt == "item":
                self.importance = clamp(self.importance+math.sqrt(self.subject.importance)/2,self.importance,self.importance*3)
                t = ""
                #if self.subject.kind in ["book","story","piece","poem","song","play"]:
                    #t = "\""
                subname = t + self.subject.name.title().capitalize().replace("\"","") + t
        if nn == "":
            n = random.choice(["a ","the ",""])
            if self.kind == "book":
                n += synonym("book")
                if random.random() > 0.3:
                    n += " " + random.choice(["about ","on ","on ","of ","relating to ","regarding "])
                else:
                    n = random.choice(["on ","about ","of ","regarding "])
                if self.subject == None:
                    if random.random() > 0.2:
                        n += synonym(self.field)
                    else:
                        n = synonym(self.field)
                else:
                    if random.random() > 0.4:
                        n += subname
                    else:
                        n += synonym(self.field)
            if self.kind in ["story","piece","poem","play","song"]:
                roll = random.random()
                if roll > 0.5:
                    if self.kind == "piece":
                        n += synonym(self.subkind)
                    else:
                        n += synonym(self.kind)
                    n += " " + random.choice(["on ","of ","regarding ","about "])
                    roll2 = random.random()
                    if roll2 > 0.4 and self.subject != None:
                        n += subname
                    elif roll2 > 0.6 and self.subject == None:
                        n += synonym(self.field)
                    else:
                        n += self.culture.language.genName()
                        if random.random() > 0.5:
                            n += " " + self.culture.language.genName()
                elif roll > 0.25:
                    n = random.choice(["the ","on ","about ","of "])
                    roll2 = random.random()
                    if roll2 > 0.2:
                        n += self.culture.language.genName()
                        if roll2 > 0.9:
                            n += " " + synonym(self.field)
                        elif roll2 > 0.75:
                            n += " " + self.culture.language.genName()
                    else:
                        n += synonym(self.field)
                        if roll2 < 0.07:
                            n += " " + self.culture.language.genName()
                else:
                    roll2 = random.random()
                    if roll2 > 0.2:
                        n = self.culture.language.genName()
                        if roll2 > 0.9:
                            n += " " + synonym(self.field)
                        elif roll2 > 0.75:
                            n += " " + self.culture.language.genName()
                    else:
                        n = synonym(self.field)
                        if roll2 < 0.07:
                            n += " " + self.culture.language.genName()
            if self.kind in ["weapon","helmet","bodice","shield","tool"]:
                roll = random.random()
                if self.kind in ["tool","weapon"]:
                    self.subkind = synonym(self.kind,exclusive=1)
                else:
                    self.subkind = synonym(self.kind)
                n = "the "
                if roll > 0.3:
                    n += self.culture.language.genName()
                elif roll > 0.6:
                    n = self.culture.language.genName()
                else:
                    n += self.subkind
                    n += " of "
                    if self.subject != None and random.random() > 0.2:
                        n += subname
                    else:
                        n += self.culture.language.genName()
        else:
            n = nn
        self.name = string.capwords(n).title()
    def move(self,loc):
        if loc == self.location:
            return -1
        if self.location != None:
            if self in self.location.items:
                self.location.items.remove(self)
        self.location = loc
        loc.items.append(self)
    def damage(self,amount,actors=[]):
        self.condition = clamp(self.condition - amount,0,1)
        if self.condition <= 0:
            e = Event(m=self.culture.myMap,a=-1,kind="destruction",sub=self,actrs=actors,loc=self.location)
            self.destructionEvent = e
            if self.location != None:
                self.location.items.remove(self)
            self.location = None
            if self.owner != None:
                self.owner.inventory.remove(self)
            self.owner = None
    def generateBookCover(self):
        self.filename = "./generated/book_"+self.name+".gif"
        sd = seedNum(self.name)
        random.seed(sd)
        marginSize = 8
        self.xDim = BookTools.bookWidth+(marginSize*2)
        self.yDim = BookTools.bookHeight+(marginSize*2)
        img = Image.new('RGB',(self.xDim,self.yDim),(255,255,255))
        backCol = (255,255,255)
        drawer = ImageDraw.Draw(img)
        numSigils = random.choice([0,0,1,1,2,2,2,3,3])
        imageCenter = (self.xDim/2,self.yDim/2)
        sigilLocation = imageCenter
        if random.random() < 0.6:
            sigilLocation = (sigilLocation[0],sigilLocation[1]-random.randint(-10,54))
        if random.random() < 0.33:
            sigilLocation = (sigilLocation[0]+random.randint(-30,30),sigilLocation[1])
        jacketColor = random.choice(BookTools.jacketColors)
        inkColor = random.choice(BookTools.inkColors)
        while inkColor == jacketColor:
            inkColor = random.choice(BookTools.inkColors)
        trimColor = random.choice(BookTools.trimColors)
        while trimColor == jacketColor:
            trimColor = random.choice(BookTools.trimColors)
        topY = marginSize
        leftX = marginSize
        bottomY = self.yDim-marginSize
        rightX = self.xDim-marginSize
        drawer.rectangle([(marginSize,marginSize),(self.xDim-marginSize,self.yDim-marginSize)],fill=jacketColor,outline=jacketColor)
        for sigilIndex in range(numSigils):
            sigilSeed = sd*(sigilIndex+1)
            newSigil = Sigil(drawer,sigilLocation,inkColor,a=-90,s=sigilSeed,o=None)
        random.seed(sd)
        namePlate = random.choice([0,1,2])
        chosenFontPath = random.choice(self.culture.language.fontPaths)
        if namePlate > 0:
            if random.random() < 0.4:
                inkColor = random.choice(BookTools.inkColors)
                while inkColor == jacketColor:
                    inkColor = random.choice(BookTools.inkColors)
            translatedTitle = self.culture.language.translatePassage(self.justName())
            fontSize = BookTools.fontSize
            chosenFont = ImageFont.truetype(chosenFontPath,fontSize)
            namePlateCenter = (self.xDim/2,self.yDim-50)
            if sigilLocation[1] >= self.yDim/2:
                namePlateCenter = (self.xDim/2,50)
            namePlateSize = drawer.multiline_textsize(translatedTitle, font=chosenFont)
            lines = 1
            while namePlateSize[0] > self.xDim-50:
                fontSize -= 1
                lines *= 2
                indexToReplace = 0
                secondIndexToReplace = 0
                while indexToReplace < len(translatedTitle)/lines:
                    indexToReplace = translatedTitle.find(" ",indexToReplace+1)
                if lines > 2:
                    secondIndexToReplace = translatedTitle.rfind(" ",0,len(translatedTitle)-indexToReplace)
                translatedTitle = translatedTitle[:indexToReplace] + "\n" + translatedTitle[indexToReplace+1:]
                if secondIndexToReplace != 0:
                    translatedTitle = translatedTitle[:secondIndexToReplace] + "\n" + translatedTitle[secondIndexToReplace+1:]
                chosenFont = ImageFont.truetype(chosenFontPath,fontSize)
                namePlateSize = drawer.multiline_textsize(translatedTitle, font=chosenFont)
            if self.culture.language.languageDirection == "right to left":
                lines = translatedTitle.split("\n")
                translatedTitle = ""
                for line in lines:
                    lineString = line[::-1]
                    translatedTitle += lineString + "\n"
                translatedTitle = translatedTitle[0:-1]
            namePlateAnchor = (namePlateCenter[0]-(namePlateSize[0]/2),namePlateCenter[1]-(namePlateSize[1]/2))
            if namePlate > 1:
                namePlateColor = trimColor
                if random.random() < 0.66:
                    namePlateColor = random.choice(BookTools.trimColors)
                namePlateTrim = None
                if random.random() < 0.4:
                    namePlateTrim = trimColor
                while namePlateColor == jacketColor or namePlateColor == inkColor:
                    namePlateColor = random.choice(BookTools.trimColors)
                drawer.rectangle([(namePlateAnchor[0]-4,namePlateAnchor[1]-4),(namePlateCenter[0]+(namePlateSize[0]/2)+4,namePlateCenter[1]+(namePlateSize[1]/2)+4)],fill=namePlateColor,outline=namePlateTrim)
            drawer.multiline_text(namePlateAnchor,translatedTitle,fill=inkColor,font=chosenFont,align="center")
        trimThickness = random.choice([0,0,2,2,2,2,2,3,4,5,6,7,8,10,12,14,16])
        if trimThickness > 0:
            topY = marginSize
            leftX = marginSize
            bottomY = self.yDim-marginSize
            rightX = self.xDim-marginSize
            trimTopY = topY+trimThickness
            trimLeftX = leftX+trimThickness
            trimBottomY = bottomY-trimThickness
            trimRightX = rightX-trimThickness
            drawer.rectangle([(leftX,topY),(rightX,trimTopY)],fill=trimColor)
            drawer.rectangle([(leftX,trimBottomY),(rightX,bottomY)],fill=trimColor)
            if random.random() > 0.5:
                if self.culture.language.languageDirection == "right to left":
                    drawer.rectangle([(leftX,topY),(trimLeftX,bottomY)],fill=trimColor)
                else:
                    drawer.rectangle([(trimRightX,topY),(rightX,bottomY)],fill=trimColor)
            else:
                drawer.rectangle([(trimRightX,topY),(rightX,bottomY)],fill=trimColor)
                drawer.rectangle([(leftX,topY),(trimLeftX,bottomY)],fill=trimColor)
        trimType = random.choice([None,None,"circle","triangle","square","circlesquare","squaretriangle","circletriangle","circlesquaretriangle"])
        if trimType != None:
            trimSize = trimThickness+random.randint(4,18)
            trimCenters = [(leftX,topY),(rightX,topY),(leftX,bottomY),(rightX,bottomY)]
            if random.random() < 0.5:
                extraTrimRoll = random.random()
                if random.random() < 0.4:
                    trimCenters = []
                if extraTrimRoll < 0.14:
                    trimCenters.extend([(leftX,self.yDim/3),(rightX,self.yDim/3)])
                elif extraTrimRoll < 0.28:
                    trimCenters.extend([(leftX,self.yDim/3),(leftX,self.yDim/1.5),(rightX,self.yDim/1.5),(rightX,self.yDim/3),(self.xDim/2,bottomY),(self.xDim/2,topY)])
                elif extraTrimRoll < 0.42:
                    trimCenters.extend([(leftX,self.yDim/3),(self.xDim/2,topY),(rightX,self.yDim/3),(self.xDim/2,bottomY)])
                elif extraTrimRoll < 0.56:
                    trimCenters.extend([(leftX,self.yDim/2),(rightX,self.yDim/2)])
                elif extraTrimRoll < 0.7:
                    trimCenters.extend([(self.xDim/2,topY),(self.xDim/2,bottomY)])
                elif extraTrimRoll < 0.84:
                    trimCenters.extend([(leftX,self.yDim/3),(rightX,self.yDim/1.5),(leftX,self.yDim/1.5),(rightX,self.yDim/3)])
                else:
                    trimCenters.extend([(leftX,self.yDim/2),(self.xDim/2,topY),(rightX,self.yDim/2),(self.xDim/2,bottomY)])
            trimIndex = 0
            for trimCenter in trimCenters:
                if trimIndex == 4:
                    if random.random() < 0.4:
                        trimType = random.choice(["circle","triangle","square"])
                        trimSize = trimThickness+random.randint(3,14)
                if "circle" in trimType:
                    drawCircle(drawer,trimCenter[0],trimCenter[1],math.floor(trimSize*1.25),trimColor) 
                if "square" in trimType:
                    drawSquare(drawer,trimCenter[0],trimCenter[1],trimSize,trimColor)
                if "triangle" in trimType:
                    drawRhombus(drawer,trimCenter[0],trimCenter[1],math.floor(trimSize*1.5),trimColor)
                trimIndex += 1
        drawer.rectangle([(0,0),(self.xDim,topY)],fill=backCol,outline=None)
        drawer.rectangle([(0,0),(leftX,self.yDim)],fill=backCol,outline=None)
        drawer.rectangle([(rightX,0),(self.xDim,self.yDim)],fill=backCol,outline=None)
        drawer.rectangle([(0,bottomY),(self.xDim,self.yDim)],fill=backCol,outline=None)
        drawer.rectangle([(marginSize,marginSize),(self.xDim-marginSize,self.yDim-marginSize)],fill=None,outline=(0,0,0))
        img.save(self.filename,"GIF")
        return self.filename
    def justName(self):
        return self.name
    def nameFull(self):
        if self.subkind != None:
            s = self.subkind
        else:
            s = self.kind
        s += " " + self.name
        if self.creator != None:
            s += " by the " + self.creator.nameFull()
        return s
    def description(self):
        vowels = Tools.vowels
        s = "\"" + self.name + "\" is a"
        if self.subkind != None:
            s = s + "n " + self.subkind if self.subkind[0].lower() in vowels else s + " " + self.subkind
        else:
            s = s + "n " + self.kind if self.kind[0].lower() in vowels else s + " " + self.kind
        if self.kind in ["weapon","helmet","bodice","shield","tool"]:
            s += " created by the " + self.creator.nameFull()
        else:
            s += " by the " + self.creator.nameFull()
        s += ".\n"
        self.material = None
        if self.kind in BookTools.writtenWorks:
            self.material = synonym("paper",seedNum(self.name))
            s += "It is written on " + self.material
            s += " in the " + self.culture.name + " language"
        elif self.subkind in ["tapestry","fresco","mural","painting","drawing"]:
            self.material = synonym("paint",seed=seedNum(self.name),exclusive=1)
            s += "It is painted with " + self.material + " paint"
        elif self.subkind in ["sculpture","statue","bust","etching"]:
            self.material = synonym("stone",seedNum(self.name))
            s += "It is made of " + self.material
        elif self.subkind in ["woodcut","longbow","shortbow","crossbow"] or self.kind in ["shield"]:
            self.material = synonym("wood",seedNum(self.name))
            s += "It is made of " + self.material
        elif self.kind in ["weapon","helmet","bodice","tool"]:
            self.material = synonym("metal",seed=seedNum(self.name))
            s += "It is made of " + self.material
        else:
            s += "It is made of mixed materials"
        if self.decoration != None:
            s += " and decorated with " + self.decoration
        s += ".\n"
        if self.kind in ["story"]:
            s += "It is fiction "
        elif self.kind == "book":
            s += "It is nonfiction "
        elif self.kind in ["piece","song","poem","play"]:
            s += "It is an art piece "
        elif self.kind == "helmet":
            s += "It is a helmet "
        elif self.kind == "bodice":
            s += "It is a piece of armor "
        elif self.kind == "shield":
            s += "It is a shield "
        elif self.kind == "weapon":
            s += "It is a weapon "
        elif self.kind == "tool":
            s += "It is a tool "
        else:
            s += "It is an object "
        if self.field != None:
            if self.kind in ["story","book","piece","poem","song","play"]:
                addition = synonym("about",seed=seedNum(self.name))
                s += addition + " " + self.field
            else:
                s += "decorated with imagery related to " + self.field
            s += ".\n"
        if self.subject != None:
            q = "work"
            if self.kind in ["helmet","weapon","bodice","shield"]:
                q = synonym("detail",seed=seedNum(self.name))
            try:
                s += "The subject of the " + q + " is the " + self.subject.culture.name + " " + self.subject.nameFull() + ".\n"
            except AttributeError: 
                s += "The subject of the " + q + " is the " + self.subject.nameFull() + ".\n"
        s += "The "
        if (self.kind in ["story","book","play","poem"]):
            s += "writing "
        elif self.kind in ["piece","song"]:
            s += "artistry "
        else:
            s += "craftsmanship "
        s += "is " + skillTier(self.quality) + ".\n"
        s += "It is generally considered "
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
        if self.condition == 0:
            s += "It is destroyed."
        return s