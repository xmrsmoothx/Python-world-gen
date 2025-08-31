# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 04:02:32 2018

@author: Bri
"""

import random
import math
import numpy as np

def factorsOf(a):
    factrs = [a]
    maximum = math.ceil(math.sqrt(a))
    for i in range(2,math.ceil(maximum)+1):
        if  a % i == 0:
            factrs.append(i)
            factrs.append(math.floor(a/i))
    factrs = list(set(factrs))
    return factrs

def drawRhombus(drawer,x,y,radius,color,out=False):
    x1 = x-radius
    x2 = x+radius
    y1 = y-radius
    y2 = y+radius
    if out == False:
        drawer.polygon([(x1,y),(x,y1),(x2,y),(x,y2)],color)
    if out == True:
        drawer.polygon([(x1,y),(x,y1),(x2,y),(x,y2)],outline=color,fill=None)

def drawSquare(drawer,x,y,radius,color,out=False):
    x1 = x-radius
    x2 = x+radius
    y1 = y-radius
    y2 = y+radius
    if out == False:
        drawer.rectangle([(x1,y1),(x2,y2)],color)
    if out == True:
        drawer.rectangle([(x1,y1),(x2,y2)],outline=color,fill=None)

def drawCircle(drawer,x,y,radius,color,out=False):
    x1 = x-radius
    x2 = x+radius
    y1 = y-radius
    y2 = y+radius
    if out == False:
        drawer.ellipse([(x1,y1),(x2,y2)],color)
    if out == True:
        drawer.ellipse([(x1,y1),(x2,y2)],outline=color,fill=None)

def lengthDirX(length, angle):
  radian_angle = math.radians(angle)
  return length * math.cos(radian_angle)

def lengthDirY(length, angle):
  radian_angle = math.radians(angle)
  return length * math.sin(radian_angle)

def getPrime(num):
    primesList = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151]
    num = num % len(primesList)
    return primesList[num]

def nearestHundred(largeNum):
    return round(largeNum/100)*100

def mode(lst):
    if lst == None or len(lst) == 0:
        return None
    lstDict = {}
    for lstElement in lst:
        if lstElement in lstDict.keys():
            lstDict[lstElement] = lstDict[lstElement]+1
        else:
            lstDict[lstElement] = 1
    mostElement = lst[0]
    for element in lstDict.keys():
        if lstDict[element] > lstDict[mostElement]:
            mostElement = element
    return mostElement

class BookTools:
    fontSize = 14
    fontPaths = ["./fonts/ManufacturingConsent-Regular.ttf",
                 "./fonts/OpenSans-Regular.ttf",
                 "./fonts/Oswald-Regular.ttf",
                 "./fonts/PTSerif-Regular.ttf",
                 "./fonts/RobotoMono-Regular.ttf",
                 "./fonts/DancingScript-Regular.ttf",
                 "./fonts/Michroma-Regular.ttf",
                 "./fonts/LobsterTwo-Regular.ttf",
                 "./fonts/Pacifico-Regular.ttf"]
    writtenWorks = ["book","story","poem","song","play"]
    bookHeight = 350
    bookWidth = 250
    jacketColors = [(70,40,40),(40,60,40),(40,40,80),(12,12,12),(60,60,60),(165,160,150),(220,190,160),(50,50,50),(46,32,21),(70,60,80),(80,70,40),(60,80,80),(235,235,235),(36,26,24)]
    inkColors = [(235,235,235),(240,180,20),(0,0,0),(46,32,21),(36,26,24)]
    trimColors = [(12,12,12),(235,235,235),(240,180,20),(60,60,60),(46,32,21),(36,26,24)]

class Tools:
    streetColor = (117, 96, 66)
    waterColor = (48, 76, 94)
    buildingColor = (40, 40, 40)
    farmColor = (156, 112, 17)
    vowels = ["a","e","i","o","u"]
    plotXDim = 600
    plotYDim = 400
    plotBackCol = (0,0,0)
    plotCol = (240,20,20)
    plotSecondCol = (0,240,0)
    mountainHeight = 0.85
    snowCapHeight = 0.91
    technologies = {}
    technologies["weaponry"] = 1
    technologies["defense"] = 1
    technologies["agriculture"] = 1
    technologies["production"] = 1
    technologies["metallurgy"] = 1
    technologies["medicine"] = 1
    technologies["government"] = 1
    technologies["art"] = 1
    technologies["philosophy"] = 1
    technologies["research"] = 1
    technologies["transportation"] = 1
    technologies["magic"] = 1

class CombatTools:
    endWarDistance = 1400
    militaryAge = 19
    warlikeSocieties = ["Pirates","Raiders"]
    aggressiveSocieties = ["Imperiun","Empire","Hegemony","Monarchy",
                                    "Religious Zealots","Revolutionary Commune",
                                    "Shamanistic Warriors"]
    battlingKinds = ["army","fleet","tactician","beast","magician"]
    armyTypes = ["assault infantry","guard infantry","siege","artillery","ranged infantry","mechanized","cavalry","fleet","beast"]
    standardArmyTypes = ["assault infantry","siege","artillery","ranged infantry","cavalry"]
    rps = {"assault infantry":["artillery","siege"],
                    "guard infantry":["assault infantry","ranged infantry","cavalry"],
                    "siege":["guard infantry","artillery"],
                    "artillery":["guard infantry","assault infantry","ranged infantry","artillery","mechanized"],
                    "ranged infantry":["assault infantry","cavalry"],
                    "mechanized":["siege","cavalry","artillery","assault infantry","ranged infantry"],
                    "cavalry":["assault infantry","artillery","siege","cavalry"],
                    }
    # How much more offensive power armies have against army types they're strong against
    rpsMultiplier = 2
    # First value is offensive modifier; second value is defensive modifier; third value is unit weight modifier
    unitBalance = {"assault infantry":[1,1,1.2],
                            "guard infantry":[1,1.5,1.2],
                            "siege":[1.25,0.75,2.5],
                            "artillery":[2,0.5,2],
                            "ranged infantry":[1.25,0.75,1.4],
                            "mechanized":[1.5,1.5,3],
                            "cavalry":[1.25,1,1.8],
                            "fleet":[10,10,4]}
    baseMilitarization = 0.08
    startingSkill = 0.4
    commandRanks = ["General","Colonel","Colonel","Commander","Commander","Commander","Captain","Captain","Captain","Captain","Captain","Captain"]
    navalRanks = ["Grand Admiral","Admiral","Admiral","Rear Admiral","Rear Admiral","Rear Admiral","Captain","Captain","Captain"]

def clamp(x,minimum,maximum):
    if x < minimum:
        return minimum
    elif x > maximum:
        return maximum
    else:
        return x

def distance3d(a,b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

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

def techTier(lv):
    steepness = 1.1
    return clamp(((2**((steepness*lv)-2))/((2**((steepness*lv)-2))+3))*(5.8*(1.002**lv))-0.4,0,8)
    

def skillTier(skill):
    tier = math.floor(skill*10)
    tiersList = ["poor","poor","fair","fair","fair","skillful","skillful","skillful","masterful","legendary","legendary"]
    return tiersList[tier]

def synonym(x,seed=0,exclusive=0):
    s = {}
    s["mountains"] = ["mountains","peaks","ridges","highlands"]
    s["savanna"] = ["savanna","plain","prairie","fields"]
    s["shrubland"] = ["shrubland","badlands","bushland"]
    s["forest"] = ["forest","woods","wood","woodland"]
    s["desert"] = ["desert","desert","wastes","barrens"]
    s["tundra"] = ["tundra","steppes","tundras"]
    s["frost tundra"] = ["frost tundra","arctic","frozen tundra","ice cap"]
    s["tropical forest"] = ["tropical forest","jungle","bush","rainforest"]
    s["boreal forest"] = ["boreal forest","woods","wood","taiga"]
    s["carnivores"] = ["carnivores","predators","hunters"]
    s["herbivores"] = ["herbivores","livestock","cattle","prey"]
    s["fear"] = ["fear","terror","dread","hatred"]
    s["warriors"] = ["warriors","fighters","soldiers"]
    s["traders"] = ["traders","merchants","sellers","brokers"]
    s["naturalism"] = ["naturalism","naturalists","herbalism","gardeners","herbalists","gardening","druids"]
    s["travelers"] = ["travelers","wanderers","nomads","migrants","itinerants","wayfarers"]
    s["sailors"] = ["sailors","mariners","wayfarers","navigation","seafaring"]
    s["swimming"] = ["swimming","sailing","navigating","cruising"]
    s["agriculture"] = ["agriculture","farming","cultivation","harvest"]
    s["camp"] = ["bivouac","camp","camp","encampment","campsite"]
    s["village"] = ["village","hamlet"]
    s["township"] = ["township","settlement"]
    s["plantlife"] = ["plantlife","plants","vegetation","flora"]
    s["vegetation"] = ["plantlife","plants","vegetation","flora"]
    s["fields"] = ["fields","farms","pastures","prairies","farmland"]
    s["metallicity"] = ["metallicity","metals","ore","prospecting","smelting"]
    s["fertility"] = ["fertility","plenty","abundance","virility","birth"]
    s["elevation"] = ["elevation","heights","mountains","cliffs"]
    s["darkness"] = ["darkness","night","twilight","dusk"]
    s["death"] = ["death","mortality","murder","the afterlife","killing"]
    s["ice"] = ["ice","snow","frost","cold"]
    s["frost"] = ["frost","permafrost","ice","glaciers","icebergs","snow"]
    s["greed"] = ["greed","wealth","gold","riches","treasure"]
    s["growth"] = ["growth","sprouting","farming","bounty","harvest"]
    s["sky"] = ["sky","stars","heavens","clouds","cosmos"]
    s["superstition"] = ["superstition","religion","faith","theology","paranormal"]
    s["collectivism"] = ["collectivism","community","cooperation","communism","socialism"]
    s["materialism"] = ["materialism","reason","empiricism","realism"]
    s["worship"] = ["worship","religion","monks","priests","prayer","mythology","clergy"]
    s["individualism"] = ["individualism","liberation","liberty","anarchism","freedom"]
    s["large"] = ["large","sizable","grand","great"]
    s["huge"] = ["huge","giant","tremendous","colossal","massive"]
    s["gigantic"] = ["gigantic","titanic","humongous","gargantuan","vast"]
    s["book"] = ["book","record","volume","document","treatise","paper","study","codex","essay","meditations","manifesto","meditation","essays","records"]
    s["story"] = ["story","novel","epic","tale","legend","chronicle"]
    s["piece"] = ["painting","woodcut","drawing","sculpture","statue","bust","etching",
     "tapestry","fresco","mural"]
    s["song"] = ["song","concerto","sonnet","ballad","opera","suite","composition","arrangement","album"]
    s["play"] = ["play","musical","opera","satire","comedy","tragedy","drama"]
    s["poem"] = ["poem","sonnet","ballad","epic"]
    s["weapon"] = ["sword","spear","greatsword","longsword","blade","rapier","crossbow","shortsword",
     "warhammer","axe","staff","sceptre","mace","lance","rifle","pistol","longbow","shortbow","halberd","pike","claymore"]
    s["helmet"] = ["helmet","helm","crown","circlet","coif","headdress","coronet","diadem","sallet","bascinet","burgonet"]
    s["bodice"] = ["bodice","breastplate","hauberk","mail","brigandine","lamellar","platemail","cuirass","coat","vest"]
    s["shield"] = ["shield","buckler","kite shield","tower shield","targe","pavise","roundshield","greatshield","small shield"]
    s["tool"] = ["tool","hammer","drill","saw","chisel","sextant","wrench","hatchet","axe","cane","brush","shovel","pickaxe"]
    s["paper"] = ["paper","parchment","vellum","slate","papyrus","bamboo","eelskin","rawhide","sandstone"]
    s["wood"] = ["wood","oak","maple","mahogany","pine","birch","hickory","fir","ash","teak","olive","cork","balsa","pecan"]
    s["stone"] = ["stone","granite","basalt","obsidian","limestone","sandstone","slate","marble","gneiss"]
    s["metal"] = ["metal","steel","iron","bronze","brass","copper","silver","gold","titanium","aluminium","tin","nickel","electrum","platinum"]
    s["paint"] = ["paint","oil","pastel","watercolor","ink","gouache","fresco","enamel","tempera"]
    s["weaponry"] = ["weaponry","combat","blades","war","battle","assault","conquest"]
    s["defense"] = ["defense","combat","armor","war","battle","fortification"]
    s["production"] = ["production","industry","manufacturing","labor"]
    s["mining"] = ["mining","minerals","mountains","metal","forging","excavation"]
    s["metallurgy"] = ["minerals","mountains","metals","forging","smithing","smelting","industry","manufacturing","labor"]
    s["government"] = ["government","bureaucracy","administration","authority","states","the state"]
    s["transportation"] = ["transportation","sailing","travel","rail","roads","infrastructure","roadbuilding"]
    s["research"] = ["research","science","experiments","physics","mathematics","language","study"]
    s["art"] = ["art","painting","sculpting","singing","music","beauty","drawing"]
    s["philosophy"] = ["philosophy","metaphysics","thought","ontology","epistemology","existentialism","knowledge","ethics"]
    s["medicine"] = ["medicine","anatomy","pharmaceuticals","surgery","illness","disease","pathogens","health"]
    s["artillery"] = ["artillery","howitzers","catapults","trebuchets","ballistas","cannons"]
    s["assault infantry"] = ["assault infantry","warriors","troopers","soldiers","infantrymen","fighters","brigade"]
    s["mechanized"] = ["mechanized","tanks","armored","engineers"]
    s["cavalry"] = ["cavalry","horseback riders","mounted","lancers","cuirassiers","horseback brigade","dragoons","hussars"]
    s["guard infantry"] = ["guard infantry","garrison","sentinels","defensive brigade","guardsmen","reserve"]
    s["ranged infantry"] = ["ranged infantry","riflemen","longbowmen","slingers","rifle brigade","carabiniers"]
    s["siege"] = ["siege","siege towers","battering rams","demolitionists","blockades","sappers"]
    s["fleet"] = ["fleet","wing","detachment","flotilla","battle group"]
    s["about"] = ["about","dealing with","related to","explaining","questioning",
     "investigating","on","concerning","relating to","challenging","exploring","pondering"]
    s["water"] = ["water","moisture","rain","rainfall","irrigation","humidity"]
    s["magic"] = ["magic","witchcraft","wizardry","miracles","sorcery","alchemy","divination","voodoo","thaumaturgy"]
    s["church"] = ["church","cathedral","parish","chapel","temple","mosque","basilica","shrine","sanctuary","abbey"]
    s["cathedral"] = ["cathedral","basilica","grand mosque","archbasilica","grand cathedral","shrine","abbey"]
    s["palace"] = ["palace","castle","citadel","keep","throne","court"]
    s["parliament"] = ["parliament","congress","capital","legislature","assembly","senate","citadel","house"]
    s["governance"] = ["governance","council","hall","assembly","town hall"]
    s["office"] = ["office","statehouse","seat","headquarters","chamber","courthouse"]
    s["longhouse"] = ["longhouse","hall","tent","meeting hall","court","fort","stronghold","grand hall"]
    s["detail"] = ["detail","engraving","filigree","embossing","brushwork","chiseling","paint"]
    s["simplicity"] = ["simplicity","asceticism","minimalism","detachment"]
    s["craftsmen"] = ["craftsmen","artisans","craft","work","creation","building","forging"]
    s["builders"] = ["builders","construction","manufacturing","craftsmen","fabrication"]
    s["ruins"] = ["ruins","wreckage"]
    s["temperature"] = ["temperature","heat","sunlight","warmth","thermodynamics"]
    s["minister"] = ["minister","secretary","chancellor","director","head","premier","chair","magistrate"]
    s["war"] = ["war","defense","security"]
    s["health"] = ["health","medicine","healthcare"]
    s["infrastructure"] = ["infrastructure","transporation"]
    s["state"] = ["state","government","the state"]
    s["culture"] = ["culture","the arts"]
    s["trade"] = ["trade","economy","exchange","finance"]
    s["party"] = ["party","league","organization","advocates","followers","guild","congress","caucus","fraternity","syndicate","conference","association","center"]
    s["cult"] = ["cult","church","clergy","worshippers","followers","priests","disciples","apostles","devotees","prophets","evangelists"]
    s["latitude"] = ["latitude","arctic","tropics"]
    s["diplomacy"] = ["diplomacy","foreign affairs","foreign policy","international relations","foreign relations"]
    s["constellations"] = ["constellations","astrology","astronomy","stars"]
    s["illness"] = ["illness","disease","flu","fever","cough","pox","sickness"]
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

class noiseField3d:
    def __init__(self,w):
        self.width = w
        self.noise = np.random.rand(w,w,w)
    def smoothNoise(self,xx,yy,zz):
        fracX = xx % 1
        fracY = yy % 1
        fracZ = zz % 1
        x1 = (math.floor(xx)+self.width) % self.width
        y1 = (math.floor(yy)+self.width) % self.width
        z1 = (math.floor(zz)+self.width) % self.width
        x2 = (x1+self.width-1) % self.width
        y2 = (y1+self.width-1) % self.width
        z2 = (z1+self.width-1) % self.width
        tileVal = 0
        tileVal += fracX*fracY*fracZ*self.noise[x1,y1,z1]
        tileVal += (1-fracX)*fracY*fracZ*self.noise[x2,y1,z1]
        tileVal += fracX*(1-fracY)*fracZ*self.noise[x1,y2,z1]
        tileVal += (1-fracX)*(1-fracY)*fracZ*self.noise[x2,y2,z1]
        tileVal += fracX*fracY*(1-fracZ)*self.noise[x1,y1,z2]
        tileVal += (1-fracX)*fracY*(1-fracZ)*self.noise[x2,y1,z2]
        tileVal += fracX*(1-fracY)*(1-fracZ)*self.noise[x1,y2,z2]
        tileVal += (1-fracX)*(1-fracY)*(1-fracZ)*self.noise[x2,y2,z2]
        return tileVal
    def turbulence(self,xx,yy,zz,size):
        tileVal = 0
        initialSize = size
        while size >= 1:
            tileVal += self.smoothNoise(xx/size,yy/size,zz/size)*size
            size = size/2
        return tileVal/(initialSize*2)

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