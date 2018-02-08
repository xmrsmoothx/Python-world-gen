# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 20:20:33 2018

@author: Bri
"""

class Event:
    def __init__(self,m=None,a=0,kind="birth",sub=None,actrs=[],loc=None):
        self.myMap = m
        self.myMap.events.append(self)
        self.age = a
        self.kind = kind
        self.subject = sub
        self.actors = actrs
        self.importance = 1
        self.location = loc
    def ageEvent(self):
        self.age += self.myMap.timeScale
    def summary(self):
        s = "the "
        s += self.kind + " of "
        s += self.subject.justName()
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
            s += " before time, "
        s += "the " + self.subject.nameFull()
        s += " was "
        if self.kind == "birth":
            s += "born"
            if self.actors == []:
                s += "."
            elif len(self.actors) == 1:
                s += " to the parent " + self.actors[0].nameFull()
            else:
                s += " to the parents "
                s += self.actors[0].nameFull() + " and "
                s += self.actors[1].nameFull() + "."
        s += "\n" + "This happened at "
        if self.location.city != None:
            c = self.location.city
            s += c.name + ", the " + c.cType(c.population)
            s += ", belonging to the " + c.culture.name + " " + c.culture.title
        elif self.location != None:
            if self.subject.culture.name not in self.location.region.culturalNames:
                s += "an unnamed " + self.location.region.biome
            else:
                s += "the "
                s += self.location.region.biome + " "
                s += self.location.region.culturalNames[self.subject.culture.name]
        s += ".\n"
        s += "This event is generally considered "
        if self.importance < 15:
            s += "unimportant"
        elif self.importance < 50:
            s += "important"
        elif self.importance < 80:
            s += "very significant"
        else:
            s += "extremely momentous"
        s += ", specifically by the "
        s += self.subject.culture.name + " culture."
        return s