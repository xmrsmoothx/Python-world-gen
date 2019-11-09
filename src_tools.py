# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 04:02:32 2018

@author: Bri
"""

import random
import math
import numpy as np

class noiseMaker:
    def __init__(self,w,h):
        self.noise = np.random.rand(w,h)
        self.width = w
        self.height = h
    def smoothNoise(self,xx,yy):
        fracX = xx % 1
        fracY = yy % 1
        x1 = (math.floor(xx)+self.width) % self.width
        y1 = (math.floor(yy)+self.height) % self.height
        x2 = (x1+self.width-1) % self.width
        y2 = (y1+self.height-1) % self.height
        tileVal = 0
        tileVal += fracX*fracY*self.noise[x1,y1]
        tileVal += (1-fracX)*fracY*self.noise[x2,y1]
        tileVal += fracX*(1-fracY)*self.noise[x1,y2]
        tileVal += (1-fracX)*(1-fracY)*self.noise[x2,y2]
        return tileVal
    def turbulence(self,xx,yy,size):
        tileVal = 0
        initialSize = size
        while size >= 1:
            tileVal += self.smoothNoise(xx/size,yy/size)*size
            size = size/2
        return tileVal/(initialSize*2)

class Tools:
    streetColor = (16,128,92)
    waterColor = (142,78,64)

def clamp(x,minimum,maximum):
    if x < minimum:
        return minimum
    elif x > maximum:
        return maximum
    else:
        return x

def discreteBlocks(x,y,s):
    return 1

def ordinal(x):
    o = str(x)
    n = x % 100
    if n == 11 or n == 12 or n == 13:
        return o + "th"
    n = x % 10
    if n == 1:
        return o + "st"
    if n == 2:
        return o + "nd"
    if n == 3:
        return o + "rd"
    else:
        return o + "th"

def synonym(x,seed=0,exclusive=0):
    s = {}
    s["mountain"] = ["mountain","peak","ridge"]
    s["savanna"] = ["savanna","plain","prairie"]
    s["shrubland"] = ["shrubland","badlands","bushland"]
    s["forest"] = ["forest","woods","wood","woodland"]
    s["desert"] = ["desert","desert","wastes","barrens"]
    s["tundra"] = ["tundra","steppes","tundra"]
    s["frost tundra"] = ["frost tundra","arctic","alpines","frozen tundra"]
    s["tropical forest"] = ["tropical forest","jungle","bush"]
    s["boreal forest"] = ["boreal forest","woods","wood","taiga"]
    s["carnivores"] = ["carnivores","predators","hunters"]
    s["herbivores"] = ["herbivores","livestock","cattle","prey"]
    s["fear"] = ["fear","terror","dread","anxiety"]
    s["warriors"] = ["warriors","fighters","soldiers"]
    s["traders"] = ["traders","merchants","sellers","brokers"]
    s["travelers"] = ["travelers","wanderers","nomads","migrants","itinerants"]
    s["agriculture"] = ["agriculture","farming","irrigation","crops","cultivation"]
    s["camp"] = ["bivouac","camp","camp","encampment","campsite"]
    s["village"] = ["village","hamlet"]
    s["township"] = ["township","settlement"]
    s["plantlife"] = ["plantlife","plants","vegetation","flora"]
    s["vegetation"] = ["plantlife","plants","vegetation","flora"]
    s["fields"] = ["fields","farms","pastures"]
    s["metallicity"] = ["metallicity","metals","ore"]
    s["fertility"] = ["fertility","plenty","abundance"]
    s["elevation"] = ["elevation","heights","mountains","cliffs"]
    s["darkness"] = ["darkness","night","twilight","dusk"]
    s["death"] = ["death","mortality","murder"]
    s["ice"] = ["ice","snow","frost","cold"]
    s["greed"] = ["greed","wealth","gold","riches"]
    s["sky"] = ["sky","stars","heavens","clouds"]
    s["superstition"] = ["superstition","religion","faith","theology","paranormal"]
    s["collectivists"] = ["collectivists","community","cooperation","communism","socialism"]
    s["freedom"] = ["freedom","liberation","liberty","anarchism"]
    s["large"] = ["large","big","sizable","grand","great"]
    s["huge"] = ["huge","giant","tremendous","colossal","massive"]
    s["gigantic"] = ["gigantic","titanic","humongous","gargantuan","vast"]
    s["book"] = ["book","record","volume","document","treatise","paper","study","codex","essay","meditations"]
    s["story"] = ["story","novel","epic","poem","tale","play","legend","chronicle"]
    s["piece"] = ["painting","woodcut","drawing","sculpture","statue","bust","etching",
     "tapestry","fresco","mural","concerto","song","sonnet","ballad"]
    s["weapon"] = ["sword","spear","greatsword","longsword","blade","rapier",
     "hammer","axe","staff","sceptre","mace","lance","rifle","pistol"]
    s["helmet"] = ["helmet","helm","crown","circlet","coif","headdress","coronet","diadem","sallet","bascinet"]
    s["bodice"] = ["bodice","breastplate","hauberk","mail coat","brigandine","lamellar","platemail","cuirass","coat"]
    s["paper"] = ["paper","parchment","vellum","slate","papyrus","bamboo"]
    s["wood"] = ["wood","oak","maple","mahogany","pine","birch","hickory","fir"]
    s["stone"] = ["stone","granite","basalt","obsidian","limestone","sandstone","slate","marble","gneiss"]
    s["metal"] = ["metal","alloy","steel","iron","bronze","brass","copper","silver","gold","titanium","aluminium","tin","nickel"]
    s["paint"] = ["paint","oil","acrylic","pastel","watercolor","ink","gouache","fresco","enamel","tempera"]
    s["weaponry"] = ["weaponry","combat","artillery","blades","war","battle","assault"]
    s["defense"] = ["defense","combat","armor","war","battle","siege","fortification"]
    s["agriculture"] = ["agriculture","farming","irrigation","crops","cultivation"]
    s["production"] = ["production","industry","factories","craftsmanship"]
    s["mining"] = ["mining","minerals","mountains","metals","forging"]
    s["metallurgy"] = ["minerals","mountains","metals","forging","smithing","smelting"]
    s["government"] = ["government","bureaucracy","administration","statism","authority","states","the state"]
    s["research"] = ["research","science","experiments","physics","mathematics","language"]
    s["equality"] = ["equality","sociology","progressivism","revolution","heirarchy","anarchy"]
    s["art"] = ["art","painting","sculpting","singing","music","beauty","drawing"]
    s["philosophy"] = ["philosophy","metaphysics","thought","ontology","epistemology","existentialism"]
    s["medicine"] = ["medicine","anatomy","pharmaceuticals","surgery","illness","disease","pathogens","health"]
    s["artillery"] = ["artillery","howitzers","catapults","trebuchets","ballistas","cannons"]
    s["assault infantry"] = ["assault infantry","warriors","troopers","soldiers","infantrymen","fighters","brigade"]
    s["mechanized"] = ["mechanized","tanks","armored","engineers"]
    s["cavalry"] = ["cavalry","horseback riders","mounted","lancers","cuirassers","horseback brigade"]
    s["guard infantry"] = ["guard infantry","garrison","sentinels","defensive brigade","guardsmen","reserve"]
    s["ranged infantry"] = ["ranged infantry","riflemen","longbowmen","slingers","rifle brigade","carbiners"]
    s["siege"] = ["siege","siege towers","battering rams","demolitionists","blockades","sappers"]
    s["about"] = ["about","dealing with","related to","explaining","questioning","investigating","on","concerning"]
    syn = x
    if x in s.keys():
        ch = random.randint(0,len(s[x])-1)
        if seed != 0:
            ch = seed % len(s[x])
        syn = s[x][ch]
    if exclusive == 1:
        if syn == x:
            syn = s[x][1]
    return syn

def seedNum(s):
    v = 0
    for k in s:
        v += ord(k)
    return v