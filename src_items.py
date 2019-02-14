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
    def __init__(self,k,c,nn="",f=None,s=None,i=1,cr=None):
        self.tt = "item"
        self.kind = k
        self.culture = c
        self.culture.items.append(self)
        self.field = f
        self.subject = s
        self.importance = i
        self.creator = cr
        self.subkind = None
        if self.kind == "piece":
            self.subkind = synonym("piece")
        subname = ""
        if self.subject != None:
            if self.subject.tt == "event":
                self.importance = clamp(self.importance+math.sqrt(self.subject.importance)/2,self.importance,self.importance*2)
                subname = self.subject.note()
            if self.subject.tt == "pop":
                self.importance = clamp(self.importance+math.sqrt(self.subject.importance)/2,self.importance,self.importance*2)
                subname = "The " + self.subject.nameFull()
            if self.subject.tt == "item":
                self.importance = clamp(self.importance+math.sqrt(self.subject.importance)/2,self.importance,self.importance*2)
                t = ""
                if self.subject.kind in ["book","story","piece"]:
                    t = "\""
                subname = t + self.subject.name + t
        if nn == "":
            n = random.choice(["a ","the "])
            if self.kind == "book":
                n += synonym("book")
                if random.random() > 0.3:
                    n += " " + random.choice(["on ","of ","relating to ","regarding "])
                else:
                    n = random.choice(["on ","about ","of ","regarding "])
                if self.subject == None:
                    if random.random() > 0.2:
                        n += synonym(self.field)
                    else:
                        n = synonym(self.field)
                else:
                    n += subname
            if self.kind == "story" or self.kind == "piece":
                roll = random.random()
                if roll > 0.65 or self.subject != None:
                    if self.kind == "story":
                        n += synonym(self.kind)
                    else:
                        n += synonym(self.subkind)
                    n += " " + random.choice(["on ","of ","regarding ","about "])
                    if self.subject == None:
                        n += self.culture.language.genName()
                        if random.random() > 0.5:
                            n += " " + self.culture.language.genName()
                    else:
                        n += subname
                elif roll > 0.3 or self.subject != None:
                    if self.subject == None:
                        if random.random() > 0.5:
                            n += synonym(self.field)
                        else:
                            n += self.culture.language.genName()
                            if random.random() > 0.5:
                                n += " " + self.culture.language.genName()
                    else:
                        n += subname
                else:
                    if self.subject == None:
                        n = self.culture.language.genName()
                        if random.random() > 0.5:
                            n += " " + self.culture.language.genName()
                    else:
                        n = subname
            if self.kind in ["weapon","helmet","bodice"]:
                roll = random.random()
                self.subkind = synonym(self.kind)
                n = "the "
                if roll > 0.3:
                    n += self.culture.language.genName()
                elif roll > 0.6:
                    n = self.culture.language.genName()
                else:
                    n += self.subkind
                    n += " of "
                    if self.subject != None:
                        n += subname
                    else:
                        n += self.culture.language.genName()
        else:
            n = nn
        self.name = string.capwords(n)
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
        s = "\"" + self.name + "\" is a "
        if self.subkind != None:
            s += self.subkind
        else:
            s += self.kind
        if self.kind == "art":
            s += " created by the " + self.creator.nameFull()
        elif self.kind in ["book","story"]:
            s += " written by the " + self.creator.nameFull()
        s += ".\n"
        self.material = None
        if (self.kind in ["story","book"] or self.subkind in 
            ["concerto","song","sonnet","ballad"]):
            self.material = synonym("paper",seedNum(self.name))
            s += "It is written on " + self.material
        elif self.subkind in ["tapestry","fresco","mural","painting","drawing"]:
            self.material = synonym("paint",seed=seedNum(self.name),exclusive=1)
            s += "It is painted with " + self.material + " paint"
        elif self.subkind in ["sculpture","statue","bust","gargoyle"]:
            self.material = synonym("stone",seedNum(self.name))
            s += "It is made of " + self.material
        elif self.subkind in ["woodcut"]:
            self.material = synonym("wood",seedNum(self.name))
            s += "It is made of " + self.material
        elif self.kind in ["weapon","helmet","bodice"]:
            self.material = synonym("alloy",seed=seedNum(self.name))
            s += "It is made of " + self.material
        else:
            s += "It is made of mixed materials"
        s += ".\n"
        if self.kind == "story":
            s += "It is fiction "
        elif self.kind == "book":
            s += "It is nonfiction "
        elif self.kind == "piece":
            s += "It is an art piece "
        else:
            s += "It is an object "
        if self.field != None:
            if self.kind in ["story","book","piece"]:
                s += "dealing with " + self.field
            else:
                s += "decorated with imagery related to " + self.field
            s += ".\n"
        if self.subject != None:
            q = "work"
            if self.kind in ["helmet","weapon","bodice"]:
                q = "engraving"
            s += "The subject of the " + q + " is the " + self.subject.nameFull() + ".\n"
        s += "It is generally considered "
        if self.importance < 12:
            s += "unimportant"
        elif self.importance < 25:
            s += "important"
        elif self.importance < 50:
            s += "extremely important"
        else:
            s += "legendary"
        s += "."
        return s