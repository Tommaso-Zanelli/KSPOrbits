#!/usr/bin/env python3
import numpy as np
import pyksp as pk
import random

# Generate new PID
def genPid(prevList = []):
    iterate = True
    while iterate:
        pid = random.randbytes(16).hex()
        if not pid in prevList:
            iterate = False
    prevList.append(pid)
    return pid

dataDict = pk.sfsFileToDict(r"C:\Users\t.zanelli\OneDrive - CINECA\quicksave #6.sfs")

pidList = []
idVessel = None
for ii, vessel in enumerate(dataDict['GAME']['FLIGHTSTATE']['VESSEL']):
    pidList.append(vessel['pid'])
    if vessel['name'] == "Around REF 0":
        idVessel = ii

vessel = dataDict['GAME']['FLIGHTSTATE']['VESSEL'][idVessel]
vessel['name'] == "Around REF 0 (Kerbol)"

# for cbId, cb in pk.orbit.system.items():
#     if cbId > 0:
#         minRad = cb.radius
#         for sat in cb.satList:
#             minRad = max(minRad, sat.orbit.apoapsis() + sat.rSOI)
#         maxRad = cb.rSOI
#         newOrb = vessel['ORBIT'].copy()
#         newOrb['SMA'] = f"{np.sqrt(minRad*maxRad):.17g}"
#         newOrb['REF'] = f"{cbId:d}"
#         newVessel = vessel.copy()
#         newVessel['pid'] = genPid(pidList)
#         newVessel['ORBIT'] = newOrb
#         newVessel['name'] = f"Around REF {cbId:d} ({cb.name:s})"
#         dataDict['GAME']['FLIGHTSTATE']['VESSEL'].append(newVessel)

#pk.dictTosfsFile(dataDict, r"C:\Users\t.zanelli\OneDrive - CINECA\testCBIDs.sfs")