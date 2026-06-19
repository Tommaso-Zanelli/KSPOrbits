#!/usr/bin/env python3
import numpy as np
import re
import json

fileId = open("rawSystemData.csv", 'rb')
csvLines = fileId.read().decode("utf-8").split('\r\n')
fileId.close()

fileIdOrder = open("rawSystemOrder.csv", 'rb')
csvLineOrders = fileIdOrder.read().decode("utf-8").split('\r\n')
fileIdOrder.close()

planets=[]
for ii in range(17):
    i1=49*ii
    i2=49*(ii + 1)
    planets.append(csvLines.copy()[i1:i2])

cbToNumDict = {}
numToCbDict = {}

for line in csvLineOrders:
    if len(line) > 0:
        lineList=line.split(',')
        refNum = int(lineList[0])
        cbName = lineList[1]
        cbToNumDict[cbName] = refNum
        numToCbDict[refNum] = cbName
    
# cbToNumDict = { "Kerbol" : 0, "Kerbin" : 1, "Mun" : 2, "Minmus" : 3, "Moho" : 4, "Eve" : 5, "Duna" : 6, "Ike" : 7, "Jool" : 8, "Laythe" : 9, "Vall" : 10, "Bop" : 11, "Tylo" : 12, "Gilly" : 13, "Pol" : 14, "Dres" : 15, "Eeloo" : 16 } 
# numToCbDict = { 0 : "Kerbol", 1 : "Kerbin", 2 : "Mun", 3 : "Minmus", 4 : "Moho", 5 : "Eve", 6 : "Duna", 7 : "Ike", 8 : "Jool", 9 : "Laythe", 10 : "Vall", 11 : "Bop", 12 : "Tylo", 13 : "Gilly", 14 : "Pol", 15 : "Dres", 16 : "Eeloo" }

def camelize(string):
    string = re.sub(r'[^\w ]', '', string)
    components = string.split()
    result = []
    for i, c in enumerate(components):
        if c.isupper():
            result.append(c)
        elif i == 0:
            result.append(c.lower())
        else:
            result.append(c.capitalize())
    return ''.join(result)

celestialBodies = {}
for ii, cb in enumerate(np.array(planets)[np.argsort([cbToNumDict[p[0]] for p in planets])]):
    cbName = cb[0]
    cbDict = {}
    cbAttList = cb[1].split(" of,")
    cbDict["class"] = cbAttList[0]
    if len(cbAttList) > 1:
        cbDict["attractor"] = { "name" : cbAttList[1], "REF" : cbToNumDict[cbAttList[1]]}
    else:
        cbDict["attractor"] = { "name" : None, "REF" : None }
    # Orbital Characteristics
    oc = cb[2:16]
    fieldName = oc[0]
    if len(fieldName) > 0:
        fieldNameCamelized = camelize(fieldName)
        ocDict = {}
        for ocfld in oc[1:]:
            oclst = ocfld.split(',')
            if len(oclst) < 1 or oclst[1] == 'None':
                ocVal = None
            else:
                ocVal = float(oclst[1])
            ocDict[camelize(oclst[0])] = ocVal
        cbDict[fieldNameCamelized] = ocDict.copy()
    else:
        cbDict['orbitalCharacteristics'] = None
    pc = cb[16:30]
    fieldName = pc[0]
    if len(fieldName) > 0:
        fieldNameCamelized = camelize(fieldName)
        pcDict = {}
        for pcfld in pc[1:]:
            pclst = pcfld.split(',')
            if len(pclst) < 1 or pclst[1] == 'None':
                pcVal = None
            else:
                pcVal = float(pclst[1])
            pcDict[camelize(pclst[0])] = pcVal
        cbDict[fieldNameCamelized] = pcDict.copy()
    ac = cb[30:41]
    fieldName = ac[0]
    if len(fieldName) > 0:
        fieldNameCamelized = camelize(fieldName)
        acDict = {}
        for acfld in ac[1:]:
            aclst = acfld.split(',')
            if len(aclst) < 1 or aclst[1] == 'None':
                acVal = None
            elif aclst[1] == 'True':
                acVal = True
            elif aclst[1] == 'False':
                acVal = False
            else:
                acVal = float(aclst[1])
            acDict[camelize(aclst[0])] = acVal
        cbDict[fieldNameCamelized] = acDict.copy() 
    sc = cb[41:]
    fieldName = sc[0]
    if len(fieldName) > 0:
        fieldNameCamelized = camelize(fieldName)
        scDict = {}
        for scfld in sc[1:]:
            sclst = scfld.split(',')
            if len(sclst) < 1 or sclst[1] == 'None':
                scVal = None
            elif sclst[1] == 'True':
                scVal = True
            elif sclst[1] == 'False':
                scVal = False
            else:
                scVal = float(sclst[1])
            scDict[camelize(sclst[0])] = scVal
        cbDict[fieldNameCamelized] = scDict.copy()
    celestialBodies[cbName] = cbDict.copy()

# Add "empty" orbit for Kerbol
celestialBodies[numToCbDict[0]]['orbitalCharacteristics'] = dict.fromkeys(celestialBodies[numToCbDict[16]]['orbitalCharacteristics'].keys(), None)

fileIdJson = open('celestialBodies.json', 'w')
json.dump(celestialBodies, fileIdJson, indent=4)
fileIdJson.close()

fileIdJson2 = open('celestialBodiesList.json', 'w')
json.dump(numToCbDict, fileIdJson2, indent=4)
fileIdJson2.close()