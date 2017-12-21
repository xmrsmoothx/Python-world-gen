import random
from Influence import Influence
from Values import Values
from Population import Population
from Flag import Flag
from Language import Language


class Culture:
    def __init__(self, n, m):
        self.origin = n
        self.myMap = m
        self.myMap.cultures.append(self)
        self.influences = Influence(self.myMap, self.origin, 1)
        self.value = Values(self.myMap, self.influences)
        self.society = self.setSociety()
        self.language = Language(self)
        self.name = self.language.name
        self.title = self.setTitle()
        self.flag = Flag(self)
        self.bannerColor = self.flag.colors[0]
        self.leaderTitle = self.setLeaderTitle()
        self.leader = Population(self, t=self.leaderTitle)

    def setSociety(self):
        m = self.value.mainValues
        if "greed" in m and "worshippers" in m and "builders" in m:
            return "Empire"
        if "warriors" in m and "greed" in m and "builders" in m and ("travelers" in m or "sailors" in m):
            return "Imperium"
        if "warriors" in m and "collectivists" in m and "worshippers" in m:
            return "Hegemony"
        if "builders" in m and "collectivists" in m and "materialists" in m:
            return "Socialists"
        if "travelers" in m and "sailors" in m and "traders" in m:
            return "Traders"
        if "freedom" in m and "collectivists" in m and "simplicity" in m:
            return "Paleolithic tribe"
        if ("travelers" in m or "sailors" in m) and ("traders" in m or "greed" in m) and "freedom" in m:
            return "Independent merchants"
        if "collectivists" in m and "agriculture" in m and "materialists" in m:
            return "Agricultural communists"
        if "collectivists" in m and "agriculture" in m and "simplicity" in m:
            return "Farming commune"
        if "worshippers" in m and "warriors" in m and ("superstitious" in m or "collectivists" in m):
            return "Religious zealots"
        if "shamans" in m and ("naturalists" in m or "astrologists" in m) and "superstitious" in m:
            return "Shamanic tribe"
        if "shamans" in m and "warriors" in m and ("astrologists" in m or "superstitious" in m or "worshippers" in m):
            return "Shamanistic warriors"
        if "metallurgists" in m and "builders" in m and "craftsmen" in m and "materialists" in m:
            return "Cooperative artisans"
        if "metallurgists" in m and "builders" in m and "craftsmen" in m and "traders" in m:
            return "Merchant artisans"
        if "freedom" in m and "greed" in m and "traders" in m:
            return "Liberal capitalists"
        if "freedom" in m and "traders" in m and ("builders" in m or "craftsmen" in m or "metallurgists" in m):
            return "Liberal merchant-artisans"
        if "builders" in m and "agriculture" in m and "traders" in m:
            return "Mercantile folk"
        if "builders" in m and "agriculture" in m and ("travelers" in m or "sailors" in m):
            return "Township builders"
        if ("craftsmen" in m or "metallurgy" in m or "builders" in m) and "simplicity" in m and (
                "naturalists" in m or "shamans" in m):
            return "Naturalist artisans"
        if ("craftsmen" in m or "metallurgy" in m or "builders" in m) and "simplicity" in m and (
                "superstitious" in m or "astrologists" in m):
            return "Traditionalist artisans"
        if "greed" in m and "sailors" in m and "warriors" in m:
            return "Pirates"
        if ("travelers" in m or "sailors" in m) and "greed" in m and "warriors" in m:
            return "Raiders"
        if ("travelers" in m or "sailors" in m) and "greed" in m and "simplicity" in m:
            return "Scavengers"
        if "travelers" in m and "simplicity" in m and "freedom" in m:
            return "Hunter-gatherer tribe"
        if "astrologists" in m and "superstitious" in m and "worshippers" in m:
            return "Religious sovereignty"
        if "collectivists" in m and "agriculture" in m and "naturalists" in m:
            return "Agriculturalists"
        if "travelers" in m and "simplicity" in m and ("metallurgists" in m or "craftsmen" in m or "builders" in m):
            return "Nomadic artisans"
        if "travelers" in m and "simplicity" and (
                "naturalists" in m or "superstitious" in m or "astrologists" in m or "shamans" in m):
            return "Nomadic tribe"
        if "materialists" in m and "metallurgists" in m and ("craftsmen" in m or "builders" in m):
            return "Blacksmiths"
        if "warriors" in m and "collectivists" in m and ("travelers" in m or "sailors" in m):
            return "Revolutionary commune"
        if "builders" in m and "metallurgists" in m and "craftsmen" in m and "agriculture" in m:
            return "Syndicalists"
        if "warriors" in m and "builders" in m and (
                "worshippers" in m or "superstitious" in m) and "freedom" not in m and "traders" not in m:
            return "Nationalists"
        if "collectivists" in m and "freedom" in m:
            return "Social Liberals"
        if (("astrologists" in m and "superstitious" in m) or
                ("superstitious" in m and "worshippers" in m) or
                ("woshippers" in m and "astrologists" in m)):
            return "Religious Sovereignty"
        if "collectivists" in m:
            return "Communalists"
        if "freedom" in m:
            return "Liberals"
        if "shamans" in m:
            return "Shamans"
        if "simplicity" in m:
            return "Tribe"
        return "Mixed society"

    def setTitle(self):
        if self.society == "Nationalists":
            return "Nation"
        if self.society == "Religious sovereignty" or self.society == "Religious zealots":
            return random.choice(["Theocracy", "Ecclesiarchy", "Order"])
        if self.society == "Agriculturalists" or self.society == "Farming commune" or self.society == "Agricultural communists":
            return random.choice(["Farmers", "Yeomen", "Peasants"])
        if self.society == "Imperium" or self.society == "Hegemony" or self.society == "Empire":
            return "Empire"
        if self.society == "Nomadic artisans" or self.society == "Nomadic peoples" or self.society == "Scavengers":
            return "Nomads"
        if (
                self.society == "Liberal capitalists" or self.society == "Liberal merchant-artisans" or self.society == "Merchant artisans" or
                self.society == "Traders" or self.society == "Independent merchants" or self.society == "Mercantile folk"):
            return "Caravans"
        if (self.society == "Blacksmiths" or self.society == "Traditionalist artisans" or
                self.society == "Naturalist artisans" or self.society == "Cooperative artisans"):
            return "Artisans"
        if (self.society == "Socialists" or self.society == "Syndicalists" or self.society == "Revolutionary commune"
                or self.society == "Communalists"):
            return random.choice(["People's Union", "Union", "Collective"])
        if (self.society == "Shamanistic warriors" or self.society == "Shamanic tribe"
                or self.society == "shamans"):
            return "Mystics"
        if self.society == "Pirates" or self.society == "Raiders":
            return "Brigands"
        if self.society == "Social Liberals" or self.society == "Liberals":
            return "Republic"
        return "People"

    def setLeaderTitle(self):
        if self.society == "Nationalists":
            return (random.choice(["Master", "High", "Lord", ""]) + " " +
                    random.choice(["Commissioner", "Chancellor", "Harbinger"]))
        if self.society == "Religious sovereignty" or self.society == "Religious zealots":
            return (random.choice(["Grand", "High", "Supreme", "Holy"]) + " " +
                    random.choice(["Pontiff", "Priest", "Shepherd"]))
        if self.society == "Agriculturalists" or self.society == "Farming commune" or self.society == "Agricultural communists":
            return (random.choice(["Head", "Chief", "Master"]) + " " +
                    random.choice(["Farmer", "Agronomist", "Foreman"]))
        if self.society == "Imperium" or self.society == "Hegemony" or self.society == "Empire":
            return "Emperor"
        if self.society == "Nomadic artisans" or self.society == "Nomadic peoples" or self.society == "Scavengers":
            return (random.choice(["Chief", "Head", "Elder"]) + " " +
                    random.choice(["Captain", "Dignitary", "Herald"]))
        if (
                self.society == "Liberal capitalists" or self.society == "Liberal merchant-artisans" or self.society == "Merchant artisans" or
                self.society == "Traders" or self.society == "Independent merchants" or self.society == "Mercantile folk"):
            return (random.choice(["Primary", "Head", "Chief", ""]) + " " +
                    random.choice(["Executive", "Director", "Superintendent"]))
        if (self.society == "Blacksmiths" or self.society == "Traditionalist artisans" or
                self.society == "Naturalist artisans" or self.society == "Cooperative artisans"):
            return (random.choice(["Master", "Elder", "Grandmaster"]) + " " +
                    random.choice(["Atificer", "Builder", "Craftsman", ""]))
        if (self.society == "Socialists" or self.society == "Syndicalists" or self.society == "Revolutionary commune"
                or self.society == "Communalists"):
            return (random.choice(["Prime", "Chief", "Central", ""]) + " " +
                    random.choice(["Director", "Governer", "Speaker"]))
        if (self.society == "Shamanistic warriors" or self.society == "Shamanic tribe"
                or self.society == "shamans"):
            return (random.choice(["Elder", "High", "Grand", "Ancestral"]) + " " +
                    random.choice(["Medicine Man", "Seer", "Shaman"]))
        if self.society == "Pirates" or self.society == "Raiders":
            return (random.choice(["Chief", "Head", ""]) + " " +
                    random.choice(["Captain", "Commander", ""]))
        if self.society == "Social Liberals" or self.society == "Liberals":
            return random.choice(["President", "Speaker", "Minister"])
        return "Chief"

    def shortName(self):
        name = ""
        name += self.name + " " + self.title
        return name

    def information(self):
        info = ""
        info += self.name + " " + self.title + "\n"
        info += "(" + self.society + ")" + "\n"
        return info

    def cultureNotes(self):
        s = self.name + " " + self.title + "\n\n"
        s += "Society type: " + self.society + "\n\n"
        s += "Leader: " + self.leader.nameFull() + "\n\n"
        s += "Capital: " + self.origin.city.name + "\n\n"
        return s
