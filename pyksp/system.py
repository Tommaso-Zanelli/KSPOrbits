#!/usr/bin/env python3
import numpy as np
from .orbit import attractor
from .orbit import orbit
from .tools import HUGE
from pathlib import Path
import json

class celestialBody(attractor):
    def __init__(self, mu_ = 1.0, rSOI_ = HUGE, name_= "", orbit_ = None, satList_ = []):
        self.mu      = mu_
        self.name    = name_
        self.orbit   = orbit_
        self.rSOI    = rSOI_
        self.satList = satList_
        if self.rSOI == HUGE:
            self.rSOI = self.computeSOI()
    
    def computeSOI(self):
        if self.orbit != None:
            mG = self.mu
            MG = self.orbit.att.mu
            a  = self.orbit.semiMajorAxis()
            return a*np.pow(mG/MG, 0.4)
        return HUGE

    def addSatellite(self, sat):
        self.satList.append(sat)

def loadSystem() -> dict[str, celestialBody]:
    with open(Path(__file__).parent.parent / "celestialBodiesList.json") as f:
        celestialBodiesList = json.load(f)
    with open(Path(__file__).parent.parent / "celestialBodies.json") as f:
        celestialBodies = json.load(f)
    system_ = {}
    for refID, cbName in celestialBodiesList.items():
        refIDNum = int(refID)
        mu_ = celestialBodies[cbName]['physicalCharacteristics']['standardGravitationalParameter']
        orbit_ = None
        attId = celestialBodies[cbName]['attractor']['REF']
        if attId != None:
            if attId in system_:
                orbit_ = orbit.jsonDictToOrbit(celestialBodies[cbName], attractor = system_[attId])
            else:
                print(f"ERROR: celestial body {cbName:s} loaded before its attractor {celestialBodies[cbName]['attractor']['name']:s}")
        system_[refIDNum] = celestialBody(mu_, HUGE, cbName, orbit_, [])
    for refIDNum, cb in system_.items():
        if cb.orbit != None:
            cb.orbit.att.satList.append(cb)
    return system_

# Module-level — loaded once on first import
SYSTEM: dict[str, celestialBody] = loadSystem()
