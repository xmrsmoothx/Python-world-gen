# -*- coding: utf-8 -*-
"""
Created on Sat Dec 25 21:04:55 2021

@author: Bri
"""

import random
from src_tools import *
from src_events import *
import string

class Magic:
    def __init__(self,c):
        self.tt = "magic"
        # Choose a type of magic spell. This is (mostly?) cosmetic.
        kinds = ["incantation","meditation","spell","prayer","invocation","channeling","concoction","ritual","song","divination"]
        # Choose one effect to do to the target; greater magnitudes are harder/less likely to be cast and generated
        effects = {"curse":-1,"bless":1,"destroy":-2.5,"create":2,"transmute":-0.25,"transport":-0.25,"damage":-0.75,"heal":0.75,"resurrect":3.5}
        # Choose a target. Greater magnitudes are harder/less likely to be cast and generated
        targets = {"item":2,"person":3,"group":4,"bloodline":5,"location":6,"city":8,"region":12,"nation":21}
        # These combinations of effects and targets won't be allowed.
        impossibleSpells = ["create region","create location",
                            "transport nation","transport location","transport region",
                            "transmute nation","transmute region","transmute location",
                            "transmute city"]
        self.kind = random.choice(kinds)
        self.effect = random.choice(list(effects.keys()))
        self.target = random.choice(list(targets.keys()))
        while self.effect + " " + self.target in impossibleSpells or abs(effects[self.effect]*targets[self.target]) > random.uniform(0,120):
            self.effect = random.choice(list(effects.keys()))
            self.target = random.choice(list(targets.keys()))
        self.strength = random.random()
        self.creator = c
        self.creator.magic.append(self)
        self.culture = self.creator.culture
        suffixes = {"curse":["doom","curse","hate","hex","spite"],
                 "bless":["blessing","sanctity","beatitude","consecration","purification"],
                 "destroy":["death","fire","inferno","mortality","horror","finality","nightmare","doom","destruction"],
                 "create":["genesis","forge","creation","primality","molding","conjuration"],
                 "transmute":["transmutation","transformation","alchemy","recreation"],
                 "transport":["teleportation","translocation","transportation","flicker"],
                 "damage":["fire","vitriol","brimstone","meteor","pain","blood"],
                 "heal":["blessing","purification","mending","healing","touch","light"],
                 "resurrect":["necromancy","revival","resurrection","unearthing","resuscitation","light"]}
        prefixes = {"curse":["doom","curse","hate","hex","spite"],
                 "bless":["holy","sacrosanct","consecrating","purifying"],
                 "destroy":["deadly","death","inferno","fire","skull","mortal","horror","final","nightmare","doom","destruction"],
                 "create":["genesis","primal","conjure"],
                 "transmute":[],
                 "transport":[],
                 "damage":["fire","vitriol","brimstone","meteor","pain","blood"],
                 "heal":["blessed","purification","mending","healing","light"],
                 "resurrect":["necromantic","soul","light","blessed","corpse","grave"]}
        s = self.culture.language.genName()
        roll = random.random()
        if roll < 0.25:
            s = s + " " + self.culture.language.genName()
        elif roll < 0.5 and len(suffixes[self.effect]) > 0:
            s = s + " " + random.choice(suffixes[self.effect])
        elif roll < 0.75 and len(prefixes[self.effect]) > 0:
            s = random.choice(prefixes[self.effect]) + " " + s
        else:
            s = s
        self.name = string.capwords(s)
        self.culture.magic.append(self)
    def justName(self):
        return self.name
    def nameFull(self):
        s = self.kind
        s += " " + self.name
        if self.creator != None:
            s += " by the " + self.creator.nameFull()
        return s
    def description(self):
        vowels = ["a","e","i","o","u"]
        s = self.name + " is a "
        if self.strength < 0.333:
            s += "weak"
        elif self.strength < 0.666:
            s += "strong"
        else:
            s += "powerful"
        s += " magic "
        s += self.kind
        s += " to "
        s += self.effect
        s += " a"
        s = s + "n " + self.target if self.target[0].lower() in vowels else s + " " + self.target
        s += ".\n"
        s += "It was created by the " + self.creator.nameFull() + "."
        return s