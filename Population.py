class Population:
    def __init__(self, c, n=None, t=""):
        self.culture = c
        if n is None:
            self.name = (self.culture.language.genName(), self.culture.language.genName())
        else:
            self.name = n
        self.title = t + " "
        self.fullName = self.nameFull()

    def nameFull(self):
        return self.title + self.name[0] + " " + self.name[1]
