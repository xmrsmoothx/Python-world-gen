from Tools import *


class Language:
    def __init__(self, c):
        self.culture = c
        self.characters()
        self.lengthPref = random.choice([3, 5, 9])
        self.name = self.genName()

    def characters(self):
        c = ['b', 'c', 'd', 'f', 'g', 'h', 'j', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'x', 'y', 'z']
        v = ['a', 'e', 'i', 'o', 'u']
        self.langConsonants = []
        self.langVowels = []
        count = 32
        n = 1
        while len(c) > 0:
            cons = random.choice(c)
            for l in range(int(math.floor(count))):
                self.langConsonants.append(cons)
            c.remove(cons)
            n += 1
            count = 32 * (1 / n) * random.uniform(0.8, 1.25)
        count = 32
        while len(v) > 0:
            vow = random.choice(v)
            for l in range(int(math.floor(count))):
                self.langVowels.append(vow)
            v.remove(vow)
            n += 1
            count = 32 * (1 / n) * random.uniform(0.8, 1.25)
        self.preferredStart = random.choice(self.langConsonants + self.langVowels)

    def genName(self):
        length = random.randint(3, 9)
        length = math.floor((length + self.lengthPref) / 2)
        n = ""
        con = 0
        vow = 0
        lastchar = '1'
        for k in range(int(length)):
            ctype = random.choice(["con", "vow"])
            if vow >= 2:
                ctype = "con"
            if con >= 2:
                ctype = "vow"
            c = lastchar
            if ctype == "con":
                c = random.choice(self.langConsonants)
                vow = 0
                con += 1
            if ctype == "vow":
                c = random.choice(self.langVowels)
                con = 0
                vow += 1
            while c == lastchar and (random.random() > 0.2):
                if ctype == "con":
                    c = random.choice(self.langConsonants)
                    vow = 0
                    con += 1
                if ctype == "vow":
                    c = random.choice(self.langVowels)
                    con = 0
                    vow += 1
            if k == 0 and random.random() < 0.35:
                c = self.preferredStart
            n += c
            lastchar = c
        return n.capitalize()
