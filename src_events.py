# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 20:20:33 2018

@author: Bri
"""

class Event:
    def __init__(self,m=None,a=0,kind="birth",sub=None,actrs=[]):
        self.myMap = m
        self.myMap.events.append(self)
        self.age = 0
        self.subject = sub
        self.actors = actrs
        self.importance = 1
    def ageEvent(self):
        self.age += self.myMap.timeScale
    def summary(self):
        s = "The "
        s += self.kind + " of "
        s += self.sub.justName()
        if self.actors != []:
            s += " by " + self.actors[0].justName()
            for k in self.actors:
                if k != self.actors[0]:
                    s += " and " + k.justName()
        s += " " + str(self.age) + " years ago"
        return s