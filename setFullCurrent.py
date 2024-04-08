# Import modules
from PyFoam.Basics.DataStructures import Vector
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from collections import namedtuple
import numpy as np
import json 
import time
import math

# reading the options file 
with open('MeshOpts') as f: 
	data = f.read() 
	
# reconstructing the data as a dictionary 
meshOpts = json.loads(data)

tStart = time.time()

try: 
    CellCenterValues = np.load('CellCenterValues.txt.npy')
    nBlocksTot = len(CellCenterValues)
    print('CellCenterValues detected and loaded in:')
    print(time.time()-tStart)
    tStart = time.time()
except:
    print('No CellCenterValues detected, determining them')
    #reading points and faces dicts
    points = np.genfromtxt('constant/polyMesh/points', skip_header=20, skip_footer=2, dtype=str)
    points = np.char.strip(np.char.strip(points, '('), ')').astype(float)

    facesStr = np.genfromtxt('constant/polyMesh/faces', skip_header=20, skip_footer=2, dtype=str,filling_values='999999',invalid_raise=False)
    faces = np.zeros([len(facesStr),5],dtype=int)
    for i in range(len(faces)):
        pts = facesStr[i][0].split("(")
        faces[i][0] = int(pts[0])
        faces[i][1] = int(pts[1])
        faces[i][2] = int(facesStr[i][1])
        faces[i][3] = int(facesStr[i][2])
        faces[i][4] = int(facesStr[i][3].strip(')'))
    print("read points and faces in:")
    print(time.time()-tStart)
    tStart = time.time()
    #reading owner and next dicts
    owner = np.genfromtxt('constant/polyMesh/owner',skip_header=21, skip_footer=2,dtype=int)
    nBlocksTot = max(owner)+1

    next = np.genfromtxt('constant/polyMesh/neighbour',skip_header=21, skip_footer=2,dtype=int)

    print("read owener and next in:")
    print(time.time()-tStart)
    tStart = time.time()
    #finding the x,y,z, of each cell center, and storing it in CellCenterValues
    nFaces  = len(faces)
    nInteralFaces = len(next)
    CellCenterValues = np.zeros([nBlocksTot,4])
    #getting the face center
    for i in range (nFaces):
        faceCenter = [0,0,0]
        #Loop over each of the 4 points
        for j in range(faces[i][0]):
            #Loop over teh 3 dimentions
            for k in range(3):
                faceCenter[k] =faceCenter[k] + points[faces[i][j+1]][k]
        faceCenter = faceCenter/faces[i][0]
        #marking that the ownwer cell has one more face counted for 
        CellCenterValues[owner[i]][3] = CellCenterValues[owner[i]][3] + 1
        #adding the x,y,z, values
        for k in range(3):
            CellCenterValues[owner[i]][k] = CellCenterValues[owner[i]][k] + faceCenter[k]
            marker = 0
            #if its an internal face, also look at the neighboor (next)
            if i < nInteralFaces:
                CellCenterValues[next[i]][k] = CellCenterValues[next[i]][k] + faceCenter[k]
        #if it was an internal face, update the neighbors count also (has to be outside of the 3D loop)
        if i < nInteralFaces:
            CellCenterValues[next[i]][3] = CellCenterValues[next[i]][3] + 1
    #averaging over the number of contributors
    for i in range(nBlocksTot):
        for k in range(3):
            CellCenterValues[i][k] = CellCenterValues[i][k] / CellCenterValues[i][3]

    print("got cell centers in:")
    print(time.time()-tStart)
    tStart = time.time()

    np.save(arr=CellCenterValues,file='CellCenterValues.txt')

    print("saved cell centers in:")
    print(time.time()-tStart)
    tStart = time.time()

currentValues = [Vector(0, 0, 0) for _ in range(nBlocksTot)]
sigmaValues = [0 for _ in range(nBlocksTot)]

wireArea = math.pi * ((meshOpts['wireRad']*1e-3)**2 - (meshOpts['wireRad']*1e-3 - meshOpts['wireShellThickness']*1e-3)**2)
J = meshOpts['I'] / wireArea
print('current density in the wire is:')
print(J)
RiseRunRatio = meshOpts['deltaH'] / (3.1415*meshOpts['coilRad']*2)
hs = np.linspace(0,meshOpts['nTurns']-1,meshOpts['nTurns'])*meshOpts['deltaH']*1e-3
sigmaValues = [0.001 for _ in range(nBlocksTot)]


#setting current in some cells
for i in range(nBlocksTot):
    #currentValues[i] = Vector(CellCenterValues[i][0],CellCenterValues[i][1],CellCenterValues[i][2])
    R2 = CellCenterValues[i][0] ** 2 + CellCenterValues[i][1] ** 2
    if R2 > ((meshOpts['coilRad']*(1e-3))-(meshOpts['wireRad']*1e-3))**2 and R2 < (meshOpts['coilRad']*1e-3+(meshOpts['wireRad']*1e-3))**2:
        hOffset = math.atan2(CellCenterValues[i][1],CellCenterValues[i][0])/2/math.pi*meshOpts['deltaH']*1e-3
        for j in range(meshOpts['nTurns']):
            if CellCenterValues[i][2] > hs[j] + hOffset - meshOpts['wireRad']*1e-3 and CellCenterValues[i][2] < hs[j] + hOffset + meshOpts['wireRad']*1e-3:
                Rmag = R2 ** (1/2)
                UnitX = CellCenterValues[i][0] / Rmag
                UnitY = CellCenterValues[i][1] / Rmag
                xLocal = CellCenterValues[i][0] - (UnitX*meshOpts['coilRad']*1e-3)
                yLocal = CellCenterValues[i][1] - (UnitY*meshOpts['coilRad']*1e-3)
                zLocal = CellCenterValues[i][2] - hs[j] - hOffset
                LocalSq = xLocal**2 + yLocal**2 + zLocal**2
                if LocalSq > (meshOpts['wireRad']*1e-3 - meshOpts['wireShellThickness']*1e-3)**2 and LocalSq < (meshOpts['wireRad']*1e-3)**2:
                    if meshOpts['method'] == 'normal' or meshOpts['method'] == 'fake' :
                        currentValues[i] = Vector(-J*UnitY*(1-RiseRunRatio),J*UnitX*(1-RiseRunRatio),J*RiseRunRatio)
                        sigmaValues[i] = 0.5
                    else:
                        sigmaValues[i] = meshOpts["sigmaWP"]
print("made the JCoil output in:")
print(time.time()-tStart)
tStart = time.time()
twoOntwo = math.sqrt(2)/2
if meshOpts['workpeice'] != "none":
    workpeiceH = [0,(meshOpts['nTurns']-1)*meshOpts["deltaH"]*(1e-3)/2,0]
    if meshOpts['workpeiceDH'] == 'auto':
        cellZ = (2 * meshOpts["extraH"] + (meshOpts["deltaH"]*(meshOpts["nTurns"])))/meshOpts["nCells"][2]
        print("cell Z is:")
        print (cellZ)
        workpeiceH[0] = (-meshOpts["extraH"]-(1.25*meshOpts["deltaH"])+(10*cellZ))
        workpeiceH[2] = (meshOpts["deltaH"]*(meshOpts["nTurns"]-1.25)+meshOpts["extraH"]-(10*cellZ))
    else:
        workpeiceH[0] = workpeiceH[1] - (meshOpts['workpeiceDH']*(1e-3))
        workpeiceH[2] = workpeiceH[1] + (meshOpts['workpeiceDH']*(1e-3))
    for i in range(nBlocksTot):
        if (meshOpts['workpeice'] == 'broken') or (meshOpts['workpeice'] == 'square'):
            if (CellCenterValues[i][2] < workpeiceH[2] and CellCenterValues[i][2] > workpeiceH[0]) or meshOpts['seperate'] == 1:
                if abs(CellCenterValues[i][0]) < 1/(meshOpts['levels'][0] + 1)*twoOntwo*(meshOpts["coilRad"]-meshOpts["wireRad"])*1e-3:
                    if abs(CellCenterValues[i][1]) < 1/(meshOpts['levels'][0] + 1)*twoOntwo*(meshOpts["coilRad"]-meshOpts["wireRad"])*1e-3:
                        sigmaValues[i] = meshOpts["sigmaWP"]
                        if meshOpts['workpeice'] == 'broken':
                            if (CellCenterValues[i][0] > meshOpts["coilRad"]*(1e-3) * 2**(1/2)/4/3): 
                                if abs(CellCenterValues[i][1]) < meshOpts["coilRad"]*(1e-3) * 2**(1/2)/ 4/3:
                                    sigmaValues[i] = 0.1
        if meshOpts['workpeice'] == 'cyl':
            R2 = CellCenterValues[i][0] ** 2 + CellCenterValues[i][1] ** 2
            if R2 < (meshOpts['workpeiceR']*1e-3)**2:
                if CellCenterValues[i][2] < workpeiceH[2] and CellCenterValues[i][2] > workpeiceH[0]:
                    if meshOpts['method'] == 'normal':
                        sigmaValues[i] = meshOpts["sigmaWP"]
                    if meshOpts['method'] == 'fake':
                        sigmaValues[i] = meshOpts["sigmaFake"]


print("made the sigma output in:")
print(time.time()-tStart)
tStart = time.time()

field_data = ParsedParameterFile('0_bkp/Jcoil')
field_data.content['internalField nonuniform List<vector>'] = currentValues
field_data.content['boundaryField'] = {'Ext' : {'type':'fixedValue',
                                                         'value':'uniform (0 0 0)'}}
field_data.writeFileAs('0/Jcoil')


print("wrote Jcoil out in:")
print(time.time()-tStart)
tStart = time.time()

field_data = ParsedParameterFile('0_bkp/sigma')
field_data['internalField nonuniform List<scalar>'] = sigmaValues

field_data.writeFileAs('0/sigma')


print("wrote sigma out in:")
print(time.time()-tStart)
tStart = time.time()
