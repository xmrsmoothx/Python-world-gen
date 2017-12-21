from Tools import random


class Values:
    def __init__(self, m, inf):
        self.myMap = m
        self.influences = inf
        self.valuesOutput = {}
        for v in self.myMap.valuesOutputs.keys():
            self.valuesOutput[v] = 0
        self.translateValues()
        self.mainValues = self.valuesMain(5)

    def maxval(self):
        maxkey = 0
        for k in self.valuesOutput.keys():
            maxkey = k
        for k in self.valuesOutput.keys():
            if self.valuesOutput[k] > self.valuesOutput[maxkey]:
                maxkey = k
        return maxkey

    def valuesMain(self, n):
        mVals = {}
        for f in range(n):
            maxkey = self.maxval()
            mVals[maxkey] = self.valuesOutput[maxkey]
            del self.valuesOutput[maxkey]
        return mVals

    def translateValues(self):
        for q in self.influences.influenceOutput.keys():
            modifier = self.influences.influenceOutput[q]
            roll = random.random()
            if roll > 0.99:
                modifier = 8
            elif roll < 0.01:
                modifier = 0.01
            for v in self.myMap.values[q].keys():
                self.valuesOutput[v] += self.myMap.values[q][v] * modifier * random.uniform(0.8, 1.25)
