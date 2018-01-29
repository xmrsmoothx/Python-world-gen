from tkinter import *
from PIL import Image, ImageDraw, ImageTk


def redraw(map_instance):
    visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
    graphDraw = ImageDraw.Draw(visualAtlas)
    if map_instance.viewmode == 0:
        for tri in map_instance.triangles:
            tri.drawReal(graphDraw,map_instance.sealevel)
        for n in map_instance.atlas:
            n.drawReal(graphDraw,map_instance.sealevel)
        for l in map_instance.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw,map_instance.xDim)
        for c in map_instance.cities:
            c.drawSelf(graphDraw)
    elif map_instance.viewmode == 1:
        for tri in map_instance.triangles:
            tri.drawTerritory(graphDraw,map_instance.sealevel)
        for l in map_instance.landmasses:
            for r in l.rivers:
                r.drawRiver(graphDraw,map_instance.xDim)
        for c in map_instance.cities:
            c.drawSelf(graphDraw)
    visualAtlas = visualAtlas.convert("RGB")
    visualAtlas.save(map_instance.mapname,"GIF")
    photo = Image.open(map_instance.mapname)
    map_instance.img = ImageTk.PhotoImage(photo)
    map_instance.lbl.configure(image = map_instance.img)
    map_instance.lbl.image = map_instance.img
    if map_instance.displayNo != None:
        map_instance.displayString.set(map_instance.nodeInfo(map_instance.displayNo))
    


def drawGraph(map_instance,gui=None):
    visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
    graphDraw = ImageDraw.Draw(visualAtlas)
    for tri in map_instance.triangles:
        tri.drawGraph(graphDraw)
    for n in map_instance.atlas:
        n.drawPoint(graphDraw,1,"red")
    visualAtlas.show()


def drawElevation(map_instance,pts=0,sea=1,gui=None):
    visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
    graphDraw = ImageDraw.Draw(visualAtlas)
    for tri in map_instance.triangles:
        tri.drawElevation(graphDraw,map_instance.sealevel*sea)
    for n in map_instance.atlas:
        n.drawElevation(graphDraw)
    for l in map_instance.landmasses:
        for r in l.rivers:
            r.drawRiver(graphDraw)
    visualAtlas.show()


def drawLandmass(map_instance,pts=0,nbrs=0,gui=None):
    visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
    graphDraw = ImageDraw.Draw(visualAtlas)
    for tri in map_instance.triangles:
        tri.drawLandmass(graphDraw,map_instance.sealevel)
    for n in map_instance.atlas:
        n.drawLandmass(graphDraw,pts,nbrs)
    for l in map_instance.landmasses:
        for r in l.rivers:
            r.drawPath(graphDraw)
    visualAtlas.show()


def drawWildlife(map_instance,gui=None):
    visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
    graphDraw = ImageDraw.Draw(visualAtlas)
    for tri in map_instance.triangles:
        tri.drawWildlife(graphDraw,map_instance.sealevel)
    for l in map_instance.landmasses:
        for r in l.rivers:
            r.drawRiver(graphDraw,map_instance.xDim)
    visualAtlas.show()


def cultureInfo(map_instance):
    if map_instance.displayNo == None:
        return -1
    if map_instance.displayNo.culture == None:
        return -1
    map_instance.displayCulture = map_instance.displayNo.culture
    if map_instance.infoGui != None:
        map_instance.infoGui.destroy()
    map_instance.infoGui = Toplevel()
    photo = Image.open(map_instance.displayCulture.flag.filename)
    map_instance.flagImg = ImageTk.PhotoImage(photo)
    map_instance.flagLbl = Label(map_instance.infoGui,image=map_instance.flagImg)
    map_instance.flagLbl.config(borderwidth=32)
    map_instance.flagLbl.photo = photo
    map_instance.flagLbl.pack()
    map_instance.cultureString = StringVar()
    map_instance.cultureString.set(map_instance.displayCulture.cultureNotes())
    cdsc = Label(map_instance.infoGui,textvariable=map_instance.cultureString)
    cdsc.pack()
    b1 = Button(map_instance.infoGui,text="Mythology Info",command=map_instance.mythologyInfo)
    b1.pack(anchor=S,side=RIGHT)
    c1 = "light goldenrod"
    b1.config(bg=c1,activebackground=c1,activeforeground=c1)


def cityInfo(map_instance):
    if map_instance.displayNo == None:
        return -1
    if map_instance.displayNo.city == None:
        return -1
    map_instance.displayCity = map_instance.displayNo.city
    if map_instance.infoGui != None:
        map_instance.infoGui.destroy()
    map_instance.infoGui = Toplevel()
    map_instance.displayCity.drawTownGen()
    photo = Image.open(map_instance.displayCity.townGen.mapName)
    map_instance.townImg = ImageTk.PhotoImage(photo)
    map_instance.townLbl = Label(map_instance.infoGui,image=map_instance.townImg)
    map_instance.townLbl.config(borderwidth=2)
    map_instance.townLbl.photo = photo
    map_instance.townLbl.pack()
    map_instance.cityString = StringVar()
    map_instance.cityString.set(map_instance.displayCity.cityNotes())
    cdsc = Label(map_instance.infoGui,textvariable=map_instance.cityString)
    cdsc.pack()
    map_instance.displayCulture = map_instance.displayNo.culture
    b1 = Button(map_instance.infoGui,text="Society Info",command=map_instance.cultureInfo)
    b1.pack(anchor=S,side=RIGHT)
    c1 = "medium aquamarine"
    b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        

def drawReal(map_instance,gui=None):
    mapDimX,mapDimY=map_instance.xDim,map_instance.yDim
    visualAtlas = Image.new("HSV",(mapDimX,mapDimY),"white")
    graphDraw = ImageDraw.Draw(visualAtlas)
    for tri in map_instance.triangles:
        tri.drawReal(graphDraw,map_instance.sealevel)
    for n in map_instance.atlas:
        n.drawReal(graphDraw,map_instance.sealevel)
    for l in map_instance.landmasses:
        for r in l.rivers:
            r.drawRiver(graphDraw,map_instance.xDim)
    for c in map_instance.cities:
        c.drawSelf(graphDraw)
    if gui == None:
        visualAtlas = visualAtlas.convert("RGB")
        visualAtlas.save("map00.png","PNG")
        visualAtlas.show()
    else:
        map_instance.gui = gui
        map_instance.displayString = StringVar()
        map_instance.displayString.set("No node selected")
        map_instance.infoScales()
        desc = Label(gui,textvariable=map_instance.displayString)
        desc.pack(side=RIGHT)
        visualAtlas = visualAtlas.convert("RGB")
        map_instance.mapname = "map_" + map_instance.cultures[0].language.genName() + ".gif"
        visualAtlas.save(map_instance.mapname,"GIF")
        photo = Image.open(map_instance.mapname)
        map_instance.img = ImageTk.PhotoImage(photo)
        map_instance.lbl = Label(gui,image=map_instance.img)
        map_instance.lbl.pack()
        map_instance.lbl.bind("<Button-1>",map_instance.displayNode)
        b0 = Button(gui,text="Next Turn",command=map_instance.nextTurn)
        b0.pack(anchor=S,side=RIGHT)
        c1 = "orange red"
        b0.config(bg=c1,activebackground=c1,activeforeground=c1)
        #b1 = Button(gui,text="Society Info",command=map_instance.cultureInfo)
        #b1.pack(anchor=S,side=RIGHT)
        c1 = "medium aquamarine"
        #b1.config(bg=c1,activebackground=c1,activeforeground=c1)
        #b3 = Button(gui,text="Settlement Info",command=map_instance.cityInfo)
        #b3.pack(anchor=S,side=RIGHT)
        c1 = "aquamarine"
        #b3.config(bg=c1,activebackground=c1,activeforeground=c1)
        b2 = Button(gui,text="Change Mode",command=map_instance.changeView)
        b2.pack(anchor=S,side=RIGHT)
        c1 = "sandy brown"
        b2.config(bg=c1,activebackground=c1,activeforeground=c1)
        gui.mainloop()
        
def drawTownGen(city_instance):
    city_instance.townGen = Town(city_instance.node,city_instance.myMap,city_instance.name)
    townImg = Image.new("HSV",(city_instance.townGen.xDim,city_instance.townGen.yDim),(0,0,255))
    graphDraw = ImageDraw.Draw(townImg)
    city_instance.townGen.drawSelf(graphDraw)
    townImg = townImg.convert("RGB")
    townImg.save(city_instance.townGen.mapName,"GIF")

def draw_flag(flag_instance):
    
    img = Image.new('HSV',(flag_instance.xDim,flag_instance.yDim),flag_instance.colors[0])
    drawer = ImageDraw.Draw(img)
    
    for shape in flag_instace.shapes:
        if shape[0]=="tri":
            drawer.polygon(*shape[1],**shape[2])
        if shape[0]=="rect":
            drawer.rectangle(*shape[1],**shape[2])
        if shape[0]=="circ":
            drawCircle(drawer,*shape[1])
            
    drawer.line(flag_instance.corners+[(0,0)],fill=(0,0,0),width=8)
    flag_instance.filename = "flag_" + flag_instance.culture.name + ".gif"
    img = img.convert('RGB')
    img.save(flag_instance.filename,"GIF")


def drawCircle(drawer,x,y,radius,color):
    x1 = x-radius
    x2 = x+radius
    y1 = y-radius
    y2 = y+radius
    drawer.ellipse([(x1,y1),(x2,y2)],color)

def drawTrapezoid(drawer,x1,y1,x2,y2,r1,r2,color):
    directAngle = A(x2-x1,y2-y1)
    pAngle = directAngle-90
    pAngle2 = pAngle-180
    p1 = (x1+lengthDirX(r1,pAngle),y1+lengthDirY(r1,pAngle))
    p2 = (x1+lengthDirX(r1,pAngle2),y1+lengthDirY(r1,pAngle2))
    p3 = (x2+lengthDirX(r2,pAngle2),y2+lengthDirY(r2,pAngle2))
    p4 = (x2+lengthDirX(r2,pAngle),y2+lengthDirY(r2,pAngle))
    drawer.polygon([p1,p2,p3,p4],color,color)
    
