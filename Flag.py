from Tools import *
from PIL import Image, ImageDraw, ImageTk


class Flag:
    def __init__(self, c):
        self.culture = c
        self.xDim = 384
        self.yDim = 192
        self.colors = []
        self.newColor()
        self.genFlag()

    def newColor(self):
        h = random.randint(0, 255)
        s = random.randint(128, 255)
        v = random.randint(0, 255)
        col = (h, s, v)
        self.colors.append(col)
        return col

    def randPt(self):
        pt = (random.randint(0, self.xDim), random.randint(0, self.yDim))
        return pt

    def center(self):
        return (self.xDim / 2, self.yDim / 2)

    def addTri(self, drawer):
        col = self.newColor()
        p0 = random.choice(self.corners)
        p1 = random.choice([(abs(self.xDim - p0[0]), p0[1]), (p0[0], abs(self.yDim - p0[1]))])
        p2 = random.choice([self.center(), random.choice(self.corners)])
        drawer.polygon([p0, p1, p2], fill=col, outline=col)

    def addRect(self, drawer):
        col = self.newColor()
        p0 = random.choice([(random.randint(0, self.xDim), random.choice([0, self.yDim, self.yDim / 2])),
                            (random.choice([0, self.xDim, self.xDim / 2]), random.randint(0, self.yDim))])
        p1 = random.choice(self.corners)
        drawer.rectangle([p0, p1], fill=col, outline=col)

    def addCirc(self, drawer):
        col = self.newColor()
        p0 = random.choice([(random.randint(0, self.xDim), random.choice([0, self.yDim, self.yDim / 2])),
                            (random.choice([0, self.xDim, self.xDim / 2]), random.randint(0, self.yDim))])
        rad = random.randint(1, self.yDim)
        drawCircle(drawer, p0[0], p0[1], rad, col)

    def genFlag(self):
        img = Image.new('HSV', (self.xDim, self.yDim), self.colors[0])
        drawer = ImageDraw.Draw(img)
        numElements = random.randint(1, 4)
        self.corners = [(0, 0), (self.xDim, 0), (0, self.yDim), (self.xDim, self.yDim)]
        for i in range(numElements):
            element = random.choice(["tri", "rect", "circ"])
            if element == "tri":
                self.addTri(drawer)
            if element == "rect":
                self.addRect(drawer)
            if element == "circ":
                self.addCirc(drawer)

        self.filename = "flag_" + self.culture.name + ".gif"
        img = img.convert('RGB')
        img.save(self.filename, "GIF")
