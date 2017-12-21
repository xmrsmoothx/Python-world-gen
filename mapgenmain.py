# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 22:47:38 2017

@author: Bri
"""

from Tools import *
from tkinter import *
from scipy.spatial import Delaunay

from Node import Node
from Map import Map
from Triangle import Triangle
from GenConfiguration import *

#----------------------------------------------------------------------#            
# Let's generate a map

atlas = [Node(-1,-1),Node(mapDimX+1,-1),Node(mapDimY+1,mapDimY+1),Node(-1,mapDimY+1)]
world = Map(atlas,numNodes,mapDimX,mapDimY)

print("Generating points...")
for x in range(numNodes-4):
    nodeX = random.random()*mapDimX
    nodeY = random.random()*mapDimY
    newNode = Node(nodeX,nodeY)
    atlas.append(newNode)

npFloatAtlas = np.zeros((numNodes,2))
for q in range(len(atlas)):
    nodeX = atlas[q].x
    nodeY = atlas[q].y
    npFloatAtlas[q] = [nodeX,nodeY]

print("Triangulating...")
triangulation = Delaunay(npFloatAtlas)

trisList = triangulation.vertices
trisVerts = triangulation.points
print("Relaxing points...")
relaxLloyd(npFloatAtlas,1)
for q in range(len(npFloatAtlas)):
    nodeX = npFloatAtlas[q,0]
    nodeY = npFloatAtlas[q,1]
    atlas[q].x = nodeX
    atlas[q].y = nodeY

triangles = []
triIndex = 0
print("Building triangles...")
while triIndex < len(trisList):
    triVertsIndices = trisList[triIndex]
    newTri = Triangle(atlas[triVertsIndices[0]],atlas[triVertsIndices[1]],atlas[triVertsIndices[2]])
    triangles.append(newTri)
    triIndex += 1

print("Assigning neighbors...")
for tri in triangles:
    if tri.verts[0].isLinked(tri.verts[1]) == 0:
        tri.verts[0].link(tri.verts[1])
    if tri.verts[1].isLinked(tri.verts[2]) == 0:
        tri.verts[1].link(tri.verts[2])
    if tri.verts[2].isLinked(tri.verts[0]) == 0:
        tri.verts[2].link(tri.verts[0])

world.triangles = triangles

print("Generating terrain...")
world.perlinElevation(6)
world.elevationAdd(-0.35)
world.addRandomShape()
world.setSeaLevel(0.4)
world.cullDots()
world.clampElevation()
world.buildAllLand()
world.buildAllWater()
world.addMajorRiver(12)
world.addMinorRiver(12)
world.cullStreams()
print("Defining biomes...")
world.setBiomes()
world.buildRegions()
world.setWildlife()
world.influences()
world.values()
world.scatterCities(16)
print("Drawing map...")
root = Tk()
world.drawReal(root)