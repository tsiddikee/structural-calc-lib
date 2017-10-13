from pythreejs import *
from IPython.display import display
from ipywidgets import HTML, Text
from traitlets import link, dlink

def renderPyThreeJsMeshes(filename, wireframeDisplay = False, nameIdentifiers = None):
    
    #create pyMeshes
    (pyMeshes, camPos) = createPyMeshes(filename, wireframeDisplay, nameIdentifiers)
    
    #create camera
    camera = PerspectiveCamera(position=camPos, fov=30, children=[DirectionalLight(color='#ffffff', position=[1, 1, 1], intensity=0.1)])
    camera.up = [0,0,1]
    
    #create scene
    pyMeshes.append(AmbientLight(color='#dddddd'))
    scene = Scene(children=pyMeshes)
    
    #render
    meshRender = Renderer(camera=camera, background='white', background_opacity=1,scene=scene, controls=[OrbitControls(controlling=camera)])
    display(meshRender)

def createPyMeshes(filename, wireframeDisplay = False, nameIdentifiers = None):
    
    #list of pyMeshes
    pyMeshes_list = []
    
    #get mesh data
    (vertices_list, faces_list, colour_list, name_list, camPos) = getMeshData(filename)
    meshCount = len(name_list)
    
    #run through meshes
    for i in range(0,meshCount):
        geometry = PlainGeometry(vertices=vertices_list[i], faces=faces_list[i])
        
        #wireframe mesh
        if(wireframeDisplay):
            #create white wireframe material
            matWire = LambertMaterial(color=rgbToHex(255,255,255))
            matWire.wireframe = True
            meshWire = Mesh(geometry=geometry, material = matWire)
            pyMeshes_list.append(meshWire)
        
        #solid mesh

        #default colour
        meshColour = rgbToHex(230,230,230)

        #mesh without highlighted elements
        if(nameIdentifiers == None):        

            #if a colour has been specified in the txt file
            if(len(colour_list) == len(name_list)):
                meshColour = rgbToHex(colour_list[i][0],colour_list[i][1],colour_list[i][2])
            
        #mesh with highlighted elements
        else:
            #test mesh name to see if the mesh shall be highlighted
            highlight = False
            if(name_list[i] in nameIdentifiers):
                highlight = True
            
            #highlighted elements in cyan
            if(highlight):
                meshColour = rgbToHex(0,170,255)
        
        matSolid = LambertMaterial(color=meshColour)
        matSolid.wireframe = False
        meshSolid = Mesh(geometry=geometry, material = matSolid)
        pyMeshes_list.append(meshSolid)
        
    return (pyMeshes_list, camPos)

def getMeshData(filename):
    
    #create empty lists to store data
    vertices_list = [] #one list per mesh
    faces_list = [] #one list per mesh
    colour_list = [] #one list per mesh
    name_list = [] #one value per mesh
    camPos = [] #One position (x,y,z) for all meshes
    
    xPos = 0.0
    yPos = 0.0
    zPos = 0.0
    
    #read file and remove whitespace \n
    with open(filename) as f:
        content = f.readlines()
        
    content = [x.strip() for x in content]   
    
    #find index of new meshes
    meshIndex = []
    for i in range(0, len(content)):
        char = content[i].split(" ")
        if(char[0] == "g"):
            meshIndex.append(i)
    meshIndex.append(len(content))
            
    #run through each mesh 
    meshCount = len(meshIndex)-1
    for i in range(0,meshCount):
        #lists
        vertices = []
        faces = []
        meshColour = []
        meshName = "Mesh"
        
        #create subset for each mesh 
        data = content[meshIndex[i]:meshIndex[i+1]]
        
        for j in range(0,len(data)):
            dataItems = data[j].split(" ")
            
            if(dataItems[0] == "v"):
                vertex = [float(x) for x in dataItems[1:]]
                vertices.append(vertex)
                
                #get outer vertex position where camera will be located
                if(vertex[0]<xPos):
                    xPos = vertex[0]
                if(vertex[1]<yPos):
                    yPos = vertex[1]
                if(vertex[2]>zPos):
                    zPos = vertex[2]
                
            elif(dataItems[0] == "c"):
                meshColour = [int(x) for x in dataItems[1:]]
                
            elif(dataItems[0] == "f"):
                face = [int(x) for x in dataItems[1:]]
                faces.append(face)
                
            elif(dataItems[0] == "g"):
                meshName = dataItems[1]
        
        #add to global lists
        vertices_list.append(vertices)
        faces_list.append(faces)
        name_list.append(meshName)
        if(len(meshColour)!=0):
            colour_list.append(meshColour)
    
    #add camera position which needs to know about all meshes
    mult = 3
    camPos.append(xPos*mult)
    camPos.append(yPos*mult)
    camPos.append(zPos*mult)
    
    #return lists
    return (vertices_list, faces_list, colour_list, name_list, camPos)

def rgbToHex(r,g,b):
    hexNumber = "#"
    
    rConvert = hex(r)[-2:]
    if(rConvert[0] == "x"):
        rConvert = rConvert[1]+rConvert[1]
    hexNumber += rConvert
    
    gConvert = hex(g)[-2:]
    if(gConvert[0] == "x"):
        gConvert = gConvert[1]+gConvert[1]
    hexNumber += gConvert
        
    bConvert = hex(b)[-2:]
    if(bConvert[0] == "x"):
        bConvert = bConvert[1]+bConvert[1]
    hexNumber += bConvert
    
    return hexNumber
