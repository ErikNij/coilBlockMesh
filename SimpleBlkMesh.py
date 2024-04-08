# Import necessary modules
from PyFoam.Basics.DataStructures import Vector
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from collections import namedtuple
import numpy as np
# importing the module 
import json
import math
#There are a couple of function her . Mist of them do what they say. Most of the time you have to input the (existing) "struct" and then something called the domain. 
#Assuming a concentric ring, the domain cotains 2 elemetns 1) the outer radis of the ring 2) The height, which is in the form of a []
#Some function may ask for some more stuff. 
def makeVerts(verts,domain):
    x = [-1,+1,+1,-1]
    y = [-1,-1,+1,+1]
    for i in range(2):
        for j in range(4):
            verts.append(Vector(x[j]*domain[0],y[j]*domain[0],domain[1][i]))

def makeVertsCenter(verts,domain):
    x = [0,+1,0,-1]
    y = [-1,0,+1,0]
    for i in range(2):
        for j in range(4):
            verts.append(Vector(x[j]*domain[0],y[j]*domain[0],domain[1][i]))

def addBlock(entry,blks,cells):
    cellsV = Vector(cells[0],cells[1],cells[2])
    blks.append('hex')
    blks.append(entry)
    blks.append(cellsV)
    blks.append('simpleGrading')
    blks.append('(1 1 1)')

def addArcPt (arcs, arcIndex,arcPoints,R):
     arcs.append('arc')
     arcs.append(arcIndex[0])
     arcs.append(arcIndex[1])
     arcs.append([R*arcPoints[0][0],R*arcPoints[0][1],arcPoints[1][0]])
     arcs.append('arc')
     arcs.append(arcIndex[0]+4)
     arcs.append(arcIndex[1]+4)
     arcs.append([R*arcPoints[0][0],R*arcPoints[0][1],arcPoints[1][1]])

def reformatentries(entry, n):
     reformed = []
     i = 0
     while i < len(entry):
          addition = ''
          for j in range(n):
               if type(entry[i])==list:
                    addition = addition  + '('
                    for k in range(len(entry[i])):
                         addition = addition + str(entry[i][k]) + ' '
                    addition = addition + ')'
               else:
                    addition = addition + str(entry[i])
               i = i + 1
               addition = addition + ' '
          reformed.append(addition)
     return reformed

def checkInput(meshOpts):
     if meshOpts['wireShellThickness'] > meshOpts['wireRad']:
          print('Shell thickness is greater than the wire!')
     if meshOpts['workpeiceR'] > meshOpts['coilRad']:
          print('Workpeice cannot fit inside the coil!')
     if meshOpts['deltaH'] < 2*meshOpts['wireRad']:
          print('Not enough space between coils in the z direction!')
# reading the options file 
with open('MeshOpts') as f: 
	data = f.read() 
	
# reconstructing the data as a dictionary 
meshOpts = json.loads(data) 

checkInput(meshOpts)
blockMeshDict = ParsedParameterFile('system/blockMeshDictBkp')
blockMeshDict.content['convertToMeters'] = 1e-3

save = 1
twoOntwo = 2**(1/2) / 2
verts = []
domains = []
arcs = []
blks = []
faces = []
#Creating the central square
if meshOpts['seperate'] == 1:
     domains.append ([1/(meshOpts['levels'][0] + 1)*twoOntwo*(meshOpts["coilRad"]-meshOpts["wireRad"]),
			[-meshOpts["workpeiceDH"]/2,meshOpts["workpeiceDH"]/2]])
if meshOpts['seperate'] == 0:
     domains.append ([1/(meshOpts['levels'][0] + 1)*twoOntwo*(meshOpts["coilRad"]-meshOpts["wireRad"]),
			[-meshOpts["extraH"]-meshOpts["deltaH"],meshOpts["deltaH"]*(meshOpts["nTurns"]-1)+meshOpts["extraH"]]])

#Creating the ring(s) between the central square and the coil
for i in range (meshOpts['levels'][0]):
     domains.append([domains[0][0]*(i+2),[-meshOpts["extraH"]-meshOpts["deltaH"],meshOpts["deltaH"]*(meshOpts["nTurns"]-1)+meshOpts["extraH"]]])



#creating the verts
for domain in domains:
     makeVerts(verts,domain)
#creating the outer domain(including the ring?)
outerDomain = []
for i in range(meshOpts['levels'][1]):
     outerDomain.append([twoOntwo*(meshOpts["coilRad"]+meshOpts["wireRad"]+(i/meshOpts['levels'][1]*meshOpts['extraWidth'])),domain[1]])
print(meshOpts['levels'][1])
for ODomain in outerDomain:
     makeVerts(verts,ODomain)
#This will help us later :)
offSet = len(verts)

if meshOpts["orginMethod"] == 'triangle':
     domainMid = domains[0]
     makeVertsCenter(verts,domainMid)
if meshOpts["orginMethod"] == 'circle':
     domainMid = domains[0]
     domainMid[0] = domainMid[0]* twoOntwo * 2  #=x * sqrt(2)
     makeVertsCenter(verts,domainMid)
if meshOpts["orginMethod"] == 'square':
     domainMid = domains[0]
     domainMid[0] = domainMid[0] * 2 #=2*x
     makeVertsCenter(verts,domainMid)
#Just need to add the orgin, and save the index
orIndex = len(verts)
verts.append(Vector(0,0,domain[1][0]))
verts.append(Vector(0,0,domain[1][1]))

blockMeshDict.content['vertices'] = verts

entries = np.linspace(start=0,stop=7,num=8,dtype=int).tolist()
entriesMid = np.linspace(start=0,stop=7,num=8,dtype=int).tolist()

#Any ring(more accurally trazoid) is the same as trazoid inisde of it, plus 8 to all the vert numbers. 
#So we just define the inner 4 trapazoids, and we will add 8. Also, we just define the bottom, as the top is the same as the bottom + 4
entries = [[8,9,1,0],
           [8,0,3,11],
           [3,2,10,11],
           [1,9,10,2]]
#ArcIndex says which of the verts are connected
arcIndex = [[0,1],[3,0],[2,3],[1,2]]
#arcPoints say in which direction (negitive y, neg x, etc The domain is also inducled. )
arcPoints = [[[0,-1],[-1,0],[0,1],[1,0]],domain[1]]
#setting the (outer) radius for each ring
R = []
for i in range(meshOpts['levels'][0]):
     R.append(domains[i+1][0]/twoOntwo)
for i in range(meshOpts['levels'][1]):
     R.append(outerDomain[i][0]/twoOntwo)
#adding the block for center? &
#adding the central square arc point
if meshOpts["orginMethod"] == "old":
     print('using old method')
     nCells = meshOpts["nCells"].copy()
     nCells[0] = nCells[0] * 1
     nCells[1] = nCells[1] * 1
     addBlock(entriesMid,blks,meshOpts['WPnCells'].copy())
     #for i in range(4):
          #addArcPt(arcs,arcIndex[i],[arcPoints[0][i],arcPoints[1]],R[0]/(meshOpts['levels'][0]+1))
     faces.append([3,2,1,0])
     faces.append([4,5,6,7])
mt = meshOpts["orginMethod"]
if (mt == 'triangle') or (mt == 'circle') or (mt =='square'):
     print('new mesh type detected')
     #entriesMid = [[0,orIndex,1,offSet],[1,orIndex,2,offSet+1],[2,orIndex,3,offSet+2],[3,orIndex,0,offSet+3]]
     entriesMid = [[0,offSet,1,orIndex],[1,offSet+1,2,orIndex],[2,offSet+2,3,orIndex],[3,offSet+3,0,orIndex]]
     addition = [4,4,4,1]
     realent = np.zeros(8,dtype=int)
     facesInt1 = []
     facesInt2 = []
     newFace = np.zeros(4,dtype=int)
     for i in range(4):
          realent[:4] = entriesMid[i]
          realent[4:] = [entriesMid[i][j] + addition[j] for j in range(4)]
          faces.append(realent[:4].tolist())
          faces[-1].reverse()
          faces.append(realent[4:].tolist())
          nCells = meshOpts["nCells"].copy()
          nCells[0] = nCells[0] / 2
          nCells[1] = nCells[1] / 2
          addBlock(realent.tolist(),blks,nCells)
          if mt == 'circle':
               for i in range(4):
                    addArcPt(arcs,arcIndex[i],[arcPoints[0][i],arcPoints[1]],R[0]/(meshOpts['levels'][0]+1))
               arcIndexCircle = [[0,offSet],[offSet,1],[1,offSet+1],[offSet+1,2],[2,offSet+2],[offSet+2,3],[3,offSet+3],[offSet+3,0]]
               smol8th = math.sin(math.radians(45/2))
               large8th = math.cos(math.radians(45/2))
               arcPointsCircle = [[[-smol8th,-large8th],[smol8th,-large8th],[large8th,-smol8th],[large8th,smol8th],[smol8th,large8th],[-smol8th,large8th],[-large8th,smol8th],[-large8th,-smol8th]],domain[1]]
     if mt == 'circle':
          for i in range(8):
               addArcPt(arcs,arcIndexCircle[i],[arcPointsCircle[0][i],arcPointsCircle[1]],R[0]/(meshOpts['levels'][0]+1))
               newFace[:2] = arcIndexCircle[i]
               upperIndex = [e + 4 for e in arcIndexCircle[i]]
               upperIndex.reverse()
               newFace[2:] = upperIndex
               #Perhaps use facesInt
               facesInt1.append(newFace.tolist())
          for i in range(4):
               newFace[:2] = arcIndex[i]
               upperIndex = [e + 4 for e in arcIndex[i]]
               upperIndex.reverse()
               newFace[2:] = upperIndex
               newFace.tolist().reverse()
               #Perhaps use facesInt
               facesInt2.append(newFace.tolist())



levels = np.linspace(start=0,num=sum(meshOpts['levels']),stop=(sum(meshOpts['levels'])-1)*8,dtype=int)
#Rlevels = [meshOpts["coilRad"]+meshOpts["wireRad"],meshOpts["coilRad"]+meshOpts["wireRad"] + meshOpts['extraWidth']]

nCells = [1,1,1]
validDims = [0,1,0,1]
#looping over the 4 Trapezoids that make up a circle
for i in range(len(entries)):
     #Looping over the different shells(levels)
     for j in range(len(levels)):
          #j = len(levels)-1
          print(j)
          #adjusting the indexs as mentioned, saing in real ent, levels is similar to offSet from the old code 
          realent = np.zeros(8,dtype=int)
          realent[:4] = [e + levels[j] for e in entries[i]]
          realent[4:] = [e + levels[j]+4 for e in entries[i]]
          if R[j] == meshOpts["coilRad"]+meshOpts["wireRad"]:
               nCells = meshOpts['nCells'].copy()
               #possibly add the refinement. However currently unused
               nCells[validDims[i]] = meshOpts['nCells'][validDims[i]]*meshOpts['refinement']
               if meshOpts['seperate'] == 1:
                    addBlock(realent.tolist(),blks,nCells)
          else:
               nCells = meshOpts['nCells'].copy()
               #add the block
          if meshOpts['orginMethod'] == 'none': 
               for k in range(8):
                    if realent[k] in [0,1,2,3]:
                         realent[k] = orIndex
                    if realent[k] in [4,5,6,7]:
                         realent[k] = orIndex+1
          if meshOpts['seperate'] == 0:
               addBlock(realent.tolist(),blks,nCells)
          addArcPt(arcs, [entries[i][arcIndex[i][0]]+levels[j],entries[i][arcIndex[i][1]]+levels[j]],[arcPoints[0][i],arcPoints[1]],R[j])
     #addArcPt(arcs, arcIndex[i],[arcPoints[0][i],arcPoints[1]],R[j])

#This is def a bit stupid, however the idea here is that the blks and arcs that we have been editing aren't quite in the right format for pyFoam, so we edit it a bit to work.
blockMeshDict.content['blocks'] = reformatentries(blks,5)
blockMeshDict.content['edges'] = reformatentries(arcs,4)

#External faces
faces.append([0+levels[-1]+8,1+levels[-1]+8,5+levels[-1]+8,4+levels[-1]+8])
faces.append([1+levels[-1]+8,2+levels[-1]+8,6+levels[-1]+8,5+levels[-1]+8])
faces.append([2+levels[-1]+8,3+levels[-1]+8,7+levels[-1]+8,6+levels[-1]+8])
faces.append([3+levels[-1]+8,0+levels[-1]+8,4+levels[-1]+8,7+levels[-1]+8])
#top and bottoms
for i in range(len(levels)):
     faces.append([0+(8*i),1+(8*i),9+(8*i),8+(8*i)])
     faces.append([1+(8*i),2+(8*i),10+(8*i),9+(8*i)])
     faces.append([2+(8*i),3+(8*i),11+(8*i),10+(8*i)])
     faces.append([3+(8*i),0+(8*i),8+(8*i),11+(8*i)])
     #AND I USE A LOOP TO DO THE TOP, LOOK AT ME DRY BB (DRY = dont repeat yourself, bb= baby)
     for i in range(4):
          if meshOpts['orginMethod'] == 'none':
               for k in range(4):
                    if faces[-4][k] in [0,1,2,3]:
                         faces[-4][k] = orIndex
                    if faces[-4][k] in [4,5,6,7]:
                         faces[-4][k] = orIndex+1
          faces.append(faces[-4])
          faces[-1].reverse()
          faces[-1] = [e + 4 if e < orIndex else e+1 for e in faces[-1]]


#Some other little (re)formatting


if meshOpts['orginMethod'] == 'circle':
     bound = np.zeros(6,dtype=type([]))
else:
     bound = np.zeros(2,dtype=type([]))

if meshOpts['orginMethod'] == 'circle':
     bound[2] = 'facesInt1'
     bound[3] = {'type':'patch',
                 'faces': []}
     for face in facesInt1:
          bound[3]['faces'].append(face)

     bound[4] = 'facesInt2'
     bound[5] = {'type':'patch',
                 'faces': []}
     for face in facesInt2:
          bound[5]['faces'].append(face)

bound[0] = 'Ext'
bound[1] = {'type':'patch',
            'faces': []}
for face in faces:
     bound[1]['faces'].append(face)

bContent = []
for i in range(len(bound)):
    bContent.append(bound[i])
#blockMeshDict.content['boundary'] = bContent

blockMeshDict.writeFileAs('system/blockMeshDict')
print("blockMeshDict created successfully.")
