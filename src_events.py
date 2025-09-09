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
        if self.subject != None and self.subject.culture != None and self.kind in ["election","reformation","war","ceasefire","battle"]:
            self.oldLeaderTitle = self.subject.culture.leaderTitle
            self.oldCultureName = self.subject.culture.shortName()
    def ageEvent(self):
        self.age += self.myMap.timeScale
    def nameFull(self):
        return self.summary()
    def justName(self):
        return self.note()
    def note(self):
        s = "The "
        if self.kind in ["war","ceasefire"]:
            s += self.kind + " with "
        if self.kind == "disbanding" and len(self.actors) > 0:
            s += "defeat of "
        else:
            s += self.kind + " of "
        if self.kind in ["reformation","war","ceasefire"]:
            s += "the " + self.oldCultureName
        elif self.kind in ["founding","disbanding"] and self.subject.tt not in ["city"]:
            s += "the " + self.subject.justName()
        else:
            s += self.subject.justName()
        return s
    def summary(self):
        s = "the "
        if self.kind == "war":
            s += self.kind + " declared on the "
        elif self.kind == "ceasefire":
            s += self.kind + " signed with the "
        elif self.kind == "battle":
            s += self.kind + " of "
        elif self.kind == "disbanding":
            if len(self.actors) == 0:
                s += self.kind + " of "
            else:
                s += " defeat of the "
        else:
            s += self.kind + " of the "
        s += self.subject.nameFull()
        if self.actors != []:
            if self.kind not in ["battle","disbanding"]:
                s += " by " + self.actors[0].justName()
                for k in self.actors:
                    if k != self.actors[0]:
                        s += " and " + k.justName()
            else:
                if self.kind == "battle":
                    s += ", fought by "
                elif self.kind == "disbanding":
                    s += " by "
                actrIndex = 1
                for k in self.actors:
                    if actrIndex < 4:
                        if k != self.actors[0]:
                            s += ", "
                        if actrIndex == 3 or k == self.actors[-1]:
                            s += "and "
                        if k.kind in ["fleet","army"]:
                            s += "the "
                        s += k.justName()
                        actrIndex += 1
        if self.year > 1:
            s += " in the year " + str(self.year)
            s += " (" + str(self.age) + " years ago)"
        elif self.subject.age < self.subject.culture.mythAge:
            s += " " + str(self.age) + " years ago "
        else:
            s += " before time"
        return s
    def fullDesc(self):
        s = ""
        if self.year > 1:
            s += "In the year " + str(self.year)
            s += " (" + str(self.age) + " years ago), "
        elif self.subject.age < self.subject.culture.mythAge:
            s += str(self.age) + " years ago, "
        else:
            s += "Before time, "
        if self.kind not in ["reformation","war","battle"]:
            s += "the " + self.subject.nameFull()
        elif self.kind == "battle":
            s += "the battle of " + self.subject.nameFull()
        else:
            s += "the " + self.oldCultureName
        if self.kind == "battle":
            s += " was fought, by "
            for k in self.actors:
                if k != self.actors[0]:
                    s += ", "
                if k == self.actors[-1]:
                    s += "and "
                s += "the "
                s += k.nameFull()
        if self.kind == "founding":
            s += " was founded"
            if len(self.actors) == 1:
                s += " by the " + self.actors[0].nameFull()
            else:
                s += " by the " + self.actors[0].nameFull()
                s += " and the " + self.actors[1].nameFull()
        if self.kind == "genesis":
            s += " came into being"
        if self.kind == "reformation":
            " was reformed into the " + self.subject.nameFull()
            if self.actors == []:
                s += ""
            else:
                s += " by the " + self.actors[0].nameFull()
        if self.kind == "war":
            s += " was attacked on the decision of the " + self.actors[0].nameFull()
        if self.kind == "ceasefire":
            s += " agreed to ceasefire with the " + self.actors[0].nameFull()
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
            s += " was elected to the position of " + self.oldLeaderTitle + " by the " + self.oldCultureName
            s += " during the election of " + str(self.year)
        if self.kind == "death":
            if len(self.actors) == 0:
                s += " died"
            else:
                s += " was killed by the " + self.actors[0].nameFull()
        if self.kind == "disbanding":
            if len(self.actors) == 0:
                s += " disbanded"
            elif len(self.actors) == 1:
                s += " was defeated by the " + self.actors[0].nameFull()
            else:
                s += " was defeated by "
                for k in self.actors:
                    if k != self.actors[0]:
                        s += ", "
                    if k == self.actors[-1]:
                        s += "and "
                    s += "the "
                    s += k.nameFull()
        if self.kind == "destruction":
            s += " was destroyed"
            if len(self.actors) > 0:
                s += " by the " + self.actors[0].nameFull()
        if self.kind == "creation":
            s += " was created"
            if self.actors == []:
                s += ""
            #elif len(self.actors) == 1:
                #s += " by the " + self.actors[0].nameFull()
        s += ".\n"
        if self.location != None:
            s += "This happened at "
            s += self.location.nodeNotes()
            s += "\n"
        s += "This event is generally considered "
        if self.importance < 10:
            s += "minor"
        elif self.importance < 21:
            s += "notable"
        elif self.importance < 34:
            s += "significant"
        elif self.importance < 50:
            s += "major"
        elif self.importance < 71:
            s += "extremely momentous"
        else:
            s += "legendary"
        s += " by the "
        s += self.subject.culture.name + " culture."
        s += "\n"
        return s