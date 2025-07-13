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
    def __init__(self,c,n=False):
        self.tt = "magic"
        self.natural = n
        self.creator = c
        self.creator.magic.append(self)
        self.culture = self.creator.culture
        # Choose a type of magic spell. This is (mostly?) cosmetic.
        kinds = ["incantation","meditation","spell","prayer","invocation","channeling","concoction","ritual","song","divination","sorcery"]
        naturalKinds = ["breath","song","roar","excretion","stare","bite"]
        # Choose one effect to do to the target; greater magnitudes are harder and less likely to be cast and generated
        effects = {"curse":-1,"bless":1,"destroy":-2.5,"create":3,"transmute":-0.4,"transport":-0.4,
                   "harm":-0.75,"heal":0.75,"resurrect":3,"burn":-1,"freeze":-1,"poison":-1.5}
        # Choose a target. Greater magnitudes are harder/less likely to be cast and generated
        targets = {"item":2,"person":3,"group":4,"bloodline":5,"location":6,"city":9,"region":15,"nation":26}
        naturalTargets = ["person","group","bloodline","location","city"]
        # These combinations of effects and targets won't be allowed.
        impossibleSpells = ["create region","create location",
                            "transport nation","transport location","transport region",
                            "transmute nation","transmute region","transmute location",
                            "transmute city"]
        self.kind = random.choice(kinds)
        if self.natural == True:
            self.kind = random.choice(naturalKinds)
        self.effect = random.choice(list(effects.keys()))
        self.target = random.choice(list(targets.keys()))
        while (self.effect + " " + self.target in impossibleSpells or abs(effects[self.effect]*targets[self.target]) > random.uniform(0,90)
        or (self.natural == True and self.target not in naturalTargets)):
            self.effect = random.choice(list(effects.keys()))
            self.target = random.choice(list(targets.keys()))
        self.strength = random.random()**2
        if self.natural == False:
            self.strength = self.strength*math.sqrt(self.culture.tech["magic"])
        self.strength = clamp((self.strength+self.creator.skill)/2,0.05,1)
        self.magnitude = effects[self.effect]*targets[self.target]
        prefixes = {"curse":["doom","curse","hate","hex","spite"],
                 "bless":["holy","sacrosanct","consecrating","purifying"],
                 "destroy":["deadly","death","inferno","fire","skull","mortal","horror","final","nightmare","doom","agony"],
                 "create":["genesis","primal","conjure","primordial"],
                 "transmute":["transmute","reshape","alchemical"],
                 "transport":["teleport","transport","translocate"],
                 "harm":["fire","vitriol","brimstone","meteor","pain","blood"],
                 "heal":["blessed","purification","mending","healing","light"],
                 "resurrect":["necromantic","soul","light","blessed","corpse","grave"],
                 "burn":["fire","brimstone","searing","combustion","inferno"],
                 "freeze":["ice","cold","arctic","biting","frost"],
                 "poison":["poison","toxic","noxious","caustic"]}
        suffixes = {"curse":["doom","curse","hate","hex","spite"],
                 "bless":["blessing","sanctity","beatitude","consecration","purification"],
                 "destroy":["death","fire","inferno","mortality","horror","finality","nightmare","doom","destruction"],
                 "create":["genesis","forge","creation","primality","molding","conjuration"],
                 "transmute":["transmutation","transformation","alchemy","recreation"],
                 "transport":["teleportation","translocation","transportation","flicker"],
                 "harm":["fire","vitriol","brimstone","meteor","pain","blood","blast"],
                 "heal":["blessing","purification","mending","healing","touch","light"],
                 "resurrect":["necromancy","revival","resurrection","unearthing","resuscitation","light"],
                 "burn":["torch","inferno","fire","brimstone","blast","pyre"],
                 "freeze":["chill","frost","icicle","blizzard","bite"],
                 "poison":["toxin","poison","infection","decay","miasma"]}
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
        if self.natural == True:
            self.name = string.capwords(self.creator.justName() + "'s " + self.kind)
        self.culture.magic.append(self)
    def cast(self,subject,caster):
        amount = self.strength
        if self.target == "item":
            self.apply(subject,amount,caster)
        if self.target == "person":
            self.apply(subject,amount,caster)
        if self.target == "group":
            self.apply(subject,amount,caster)
        if self.target == "bloodline":
            self.apply(subject,amount,caster)
            for k in subject.getDescendants():
                self.apply(k,amount,caster)
        if self.target == "location":
            self.apply(subject,amount,caster)
        if self.target == "city":
            self.apply(subject.node,amount,caster)
        if self.target == "region":
            for n in subject.nodes:
                self.apply(n,amount,caster)
        if self.target == "nation":
            for c in subject.cities:
                rr = c.node.resourceRegion
                for n in rr.nodes:
                    self.apply(n,amount,caster)
    def apply(self,subject,amount,caster):
        if subject.tt == "item":
            if self.effect == "destroy":
                subject.damage(1,[caster])
        if subject.tt == "pop":
            if self.effect == "destroy":
                subject.kill(1,[caster])
        if subject.tt == "node":
            if self.effect == "destroy":
                a = 1
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
        s = self.name + " is a"
        if self.strength < 0.25:
            s += ""
        elif self.strength < 0.5:
            s += " strong"
        elif self.strength < 0.75:
            s += " powerful"
        else:
            s += " legendary"
        s += " magic "
        s += self.kind
        s += " to "
        s += self.effect
        s += " a"
        s = s + "n " + self.target if self.target[0].lower() in vowels else s + " " + self.target
        s += ".\n"
        if self.natural == False:
            s += "It was created by the " + self.creator.nameFull() + "."
        else:
            s += "It is an ability of the " + self.creator.nameFull() + "."
        return s