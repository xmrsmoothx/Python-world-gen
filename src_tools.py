# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 04:02:32 2018

@author: Bri
"""

import random
import math

def synonym(x,seed=0):
    s = {}
    s["mountain"] = ["mountain","peak","ridge"]
    s["savanna"] = ["savanna","plain","prairie"]
    s["shrubland"] = ["shrubland","badlands","bushland"]
    s["forest"] = ["forest","woods","wood","woodland"]
    s["desert"] = ["desert","desert","wastes","barrens"]
    s["tundra"] = ["tundra","steppes","tundra"]
    s["frost tundra"] = ["frost tundra","arctic","alpines","frozen tundra"]
    s["tropical forest"] = ["tropical forest","jungle"]
    s["boreal forest"] = ["boreal forest","woods","wood","taiga"]
    s["carnivores"] = ["carnivores","predators"]
    s["herbivores"] = ["herbivores","livestock","cattle"]
    s["fear"] = ["fear","terror"]
    s["warriors"] = ["warriors","fighters","soldiers"]
    s["agriculture"] = ["agriculture","farming","irrigation","crops","cultivation"]
    s["camp"] = ["bivouac","camp","camp","encampment","campsite"]
    s["village"] = ["village","hamlet"]
    s["township"] = ["township","settlement"]
    s["plantlife"] = ["plantlife","plants","vegetation","flora"]
    s["vegetation"] = ["plantlife","plants","vegetation","flora"]
    s["fields"] = ["fields","farms","pastures"]
    s["metallicity"] = ["metallicity","metals","ore"]
    s["fertility"] = ["fertility","plenty","abundance"]
    s["darkness"] = ["darkness","night","twilight","dusk"]
    s["death"] = ["death","mortality"]
    s["ice"] = ["ice","snow","frost"]
    s["greed"] = ["greed","wealth","gold"]
    s["sky"] = ["sky","stars","heavens"]
    s["superstition"] = ["superstition","religion","faith","theology","paranormal"]
    s["collectivists"] = ["collectivists","community","cooperation"]
    s["freedom"] = ["freedom","liberation","liberty"]
    s["large"] = ["large","big","sizable","oversized","bulky"]
    s["huge"] = ["huge","giant","enormous","colossal","immense"]
    s["gigantic"] = ["gigantic","tremendous","titanic","humongous","gargantuan"]
    s["book"] = ["book","record","volume","document","treatise","paper","study","codex","essay"]
    s["story"] = ["story","novel","epic","poem","tale","play","legend","chronicle"]
    s["piece"] = ["painting","woodcut","drawing","sculpture","statue","bust","gargoyle",
     "tapestry","fresco","mural","concerto","song","sonnet","ballad","painting"]
    s["weapon"] = ["sword","spear","greatsword","longsword","blade",
     "hammer","axe","staff","sceptre","mace","lance","rifle","pistol"]
    s["helmet"] = ["helmet","helm","crown","circlet","coif","headdress","coronet","diadem","sallet","bascinet"]
    s["bodice"] = ["bodice","breastplate","hauberk","mail coat","brigandine","lamellar","platemail"]
    s["paper"] = ["paper","parchment","vellum","slate","papyrus","bamboo"]
    s["wood"] = ["wood","oak","maple","mahogany","pine","birch","hickory","fir"]
    s["stone"] = ["stone","granite","basalt","obsidian","limestone","sandstone","slate","marble","gneiss"]
    s["alloy"] = ["alloy","steel","iron","bronze","brass","copper","silver","gold","titanium","aluminium"]
    s["paint"] = ["oil","acrylic","pastel","watercolor","ink","gouache","fresco","enamel","tempera"]
    s["weaponry"] = ["weaponry","combat","artillery","blades","war","battle"]
    s["defense"] = ["defense","combat","armor","war","battle","siege"]
    s["agriculture"] = ["agriculture","farming","irrigation","crops","cultivation"]
    s["production"] = ["production","industry","factories","craftsmanship"]
    s["mining"] = ["mining","minerals","mountains","metals","forging"]
    s["metallurgy"] = ["minerals","mountains","metals","forging","smithing"]
    s["government"] = ["government","bureaucracy","administration","statism","authority","states"]
    s["research"] = ["research","science","experiments","physics","mathematics"]
    s["equality"] = ["equality","sociology","progressivism","revolution","heirarchy","anarchy"]
    s["art"] = ["art","painting","sculpting","singing","music","beauty","drawing"]
    s["philosophy"] = ["philosophy","metaphysics","thought","ontology","epistemology","existentialism"]
    s["medicine"] = ["medicine","anatomy","pharmaceuticals","surgery","illness","disease","pathogens"]
    syn = x
    if x in s.keys():
        ch = random.randint(0,len(s[x])-1)
        if seed != 0:
            ch = seed % len(s[x])
        syn = s[x][ch]
    return syn

def seedNum(s):
    v = 0
    for k in s:
        v += ord(k)
    return v