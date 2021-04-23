# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 20:20:33 2018

@author: Bri
"""

import random
import math

class Event:
    def __init__(self,m=None,a=1,kind="birth",sub=None,actrs=[],loc=None):
        self.tt = "event"
        self.myMap = m
        self.myMap.events.append(self)
        self.age = a
        self.year = self.myMap.date-self.age
        self.kind = kind
        self.subject = sub
        self.actors = actrs
        self.importance = math.ceil(50*(random.uniform(0.1,0.6)**2))
        self.location = loc
    def ageEvent(self):
        self.age += self.myMap.timeScale
    def nameFull(self):
        return self.summary()
    def justName(self):
        return self.note()
    def note(self):
        s = "The "
        s += self.kind + " of "
        s += self.subject.justName()
        return s
    def summary(self):
        s = "the "
        s += self.kind + " of the "
        s += self.subject.nameFull()
        if self.actors != []:
            s += " by " + self.actors[0].justName()
            for k in self.actors:
                if k != self.actors[0]:
                    s += " and " + k.justName()
        if self.subject.age < self.subject.culture.mythAge:
            s += " " + str(self.age) + " years ago"
        else:
            s += " before time"
        return s
    def fullDesc(self):
        s = ""
        if self.subject.age < self.subject.culture.mythAge:
            s += "" + str(self.age) + " years ago, "
        else:
            s += "Before time, "
        s += "the " + self.subject.nameFull()
        if self.kind == "founding":
            s += " was founded"
            if len(self.actors) == 1:
                s += " by the " + self.actors[0].nameFull()
        if self.kind == "genesis":
            s += " came into being"
        if self.kind == "birth":
            s += " was born"
            if self.actors == []:
                s += ""
            elif len(self.actors) == 1:
                s += " to "+self.subject.possessive[self.subject.gender]+" parent the " + self.actors[0].nameFull()
            else:
                s += " to "+self.subject.possessive[self.subject.gender]+" parents "
                s += "the " + self.actors[0].nameFull() + " and "
                s += "the " + self.actors[1].nameFull()
        if self.kind == "election":
            s += " was elected to the leadership of the " + self.subject.culture.shortName()
            s += " during the election of " + str(self.year)
        if self.kind == "death":
            s += " died"
        if self.kind == "disbanded":
            s += " disbanded"
        s += ".\n"
        if self.location != None:
            s += "This happened at "
            if self.location.city != None:
                c = self.location.city
                s += c.name + ", the " + c.cType(c.population)
                s += ", belonging to the " + c.culture.name + " " + c.culture.title+".\n"
            else:
                if self.subject.culture.name not in self.location.region.culturalNames:
                    s += "an unnamed " + self.location.region.biome
                else:
                    s += "the "
                    s += self.location.region.biome + " "
                    s += self.location.region.culturalNames[self.subject.culture.name]
                s += ".\n"
        s += "This event is generally considered "
        if self.importance < 15:
            s += "minor"
        elif self.importance < 50:
            s += "major"
        elif self.importance < 80:
            s += "extremely momentous"
        else:
            s += "legendary"
        s += " by the "
        s += self.subject.culture.name + " culture."
        s += "\n"
        return s