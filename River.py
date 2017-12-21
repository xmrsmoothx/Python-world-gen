from Tools import *


class River:
    def __init__(self, root, length, landms):
        self.nodes = []
        self.addNode(root)
        self.landmass = landms
        current = root
        j = 0
        while j < length:
            choice = current
            for m in current.neighbors:
                if (m.elevation > choice.elevation
                        and m.river is None and m not in self.nodes):
                    choice = m
            for k in current.neighbors:
                if (k.elevation < choice.elevation and k.elevation >= current.elevation
                        and k.river is None and k not in self.nodes):
                    choice = k
            if current == choice:
                j = length
            else:
                self.addNode(choice)
            j += 1
            current = choice
        self.culturalNames = {}

    def removeRiver(self):
        for n in self.nodes:
            n.river = None
            self.nodes.remove(n)

    def addNode(self, newNode):
        self.nodes.append(newNode)
        newNode.river = self

    def drawPath(self, drawer):
        dCol = (142, 64, 64)
        drawNodes = []
        for n in self.nodes:
            drawNodes.append(n.coords())
        drawer.line(drawNodes, dCol, 2)

    def drawRiver(self, drawer, xDim):
        for l in self.nodes:
            l.getSlope()
        dCol = (142, 128, 76)
        for i in range(len(self.nodes) - 1):
            n = self.nodes[i]
            n1 = self.nodes[i + 1]
            scale = xDim / 2
            w = clamp(1 / n.slope, 0, 2048) / scale
            w1 = clamp(1 / n1.slope, 0, 2048) / scale
            drawCircle(drawer, n.x, n.y, w, dCol)
            drawCircle(drawer, n1.x, n1.y, w1, dCol)
            drawTrapezoid(drawer, n.x, n.y, n1.x, n1.y, w, w1, dCol)
