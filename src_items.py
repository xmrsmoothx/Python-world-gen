# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 01:07:42 2018

@author: Bri
"""

import random
from src_tools import *
from src_events import *
import string
import math

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
            n = random.choice(["a ","the "])
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
        if self.kind in ["story","book","poem","play","song"]:
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