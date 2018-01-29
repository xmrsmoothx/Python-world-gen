
import src_mapgenmain as src_map
import draw
from tkinter import *

import numpy as np
import random

import json

#----------------------------------------------------------------------#            
# Let's generate a map

def generate(save=True):
    
    numNodes = 2**6#2**14
    mapDimX = 960
    mapDimY = 960
    if save:
        with open("basic_info.txt","w") as f:
            d={"numNodes":numNodes,
                "mapDimX":mapDimX,
                "mapDimY":mapDimY,}
            
            d_s=json.dumps(d)
            
            f.write(d_s)
        
    atlas = [src_map.Node(-1,-1),src_map.Node(mapDimX+1,-1),src_map.Node(mapDimY+1,mapDimY+1),src_map.Node(-1,mapDimY+1)]
    world = src_map.Map(atlas,numNodes,mapDimX,mapDimY)

    print("Generating points...")
    for x in range(numNodes-4):
        nodeX = random.random()*mapDimX
        nodeY = random.random()*mapDimY
        newNode = src_map.Node(nodeX,nodeY)
        atlas.append(newNode)

    npFloatAtlas = np.zeros((numNodes,2))
    for q in range(len(atlas)):
        nodeX = atlas[q].x
        nodeY = atlas[q].y
        npFloatAtlas[q] = [nodeX,nodeY]

    print("Triangulating...")
    from scipy.spatial import Delaunay
    triangulation = Delaunay(npFloatAtlas)

    trisList = triangulation.vertices
    trisVerts = triangulation.points
    print("Relaxing points...")
    src_map.relaxLloyd(npFloatAtlas,1)
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
        newTri = src_map.Triangle(atlas[triVertsIndices[0]],atlas[triVertsIndices[1]],atlas[triVertsIndices[2]])
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
    world.soilProperties()
    print("Defining biomes...")
    world.setBiomes()
    world.buildRegions()
    world.setWildlife()
    world.influences()
    world.values()
    world.godSpheres()
    world.scatterCities(16)
    print("Drawing map...")
    root = Tk()
    draw.drawReal(world,root)
    #world.drawReal(root)
    
    
if __name__=="__main__":
    generate()
