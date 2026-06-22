#!/usr/bin/env python3
import numpy as np
import oct2py as op
from oct2py import octave as oc
from .tools import tiny, HUGE, degToRad, radToDeg
from pathlib import Path
"""
Provides an interface between python and the orbital mechanics library written    
by me for matlab a few years ago

Created on Tue Nov 12 15:18:07 2024

@author: Tommaso Zanelli
"""

path = str(Path(__file__).parent.parent / "OctaveScripts")
oc.addpath(path)
oc.addpath(oc.genpath(path))

class attractor:
    def __init__(self, mu_ = 1.0, rSOI_ = HUGE):
        self.mu   = mu_
        self.rSOI = rSOI_
    
defaultAttractor = attractor()
    
class orbit:
    system: dict[str, attractor] = {}
    
    def __init__(self):
        self.p    = 1.0
        self.e    = 0.0
        self.th0  = 0.0
        self.incl = 0.0
        self.lan  = 0.0
        self.t0   = 0.0
        self.att  = defaultAttractor

    @classmethod
    def defineParameters(cls, p_, e_ = 0.0, th0_ = 0.0, incl_ = 0.0, lan_ = 0.0, t0_ = 0.0, att_ = defaultAttractor):
        orbit_ = cls()
        orbit_.p    = p_
        orbit_.e    = e_
        orbit_.th0  = th0_
        orbit_.incl = incl_
        orbit_.lan  = lan_
        orbit_.t0   = t0_
        orbit_.att  = att_
        return orbit_

    @classmethod
    def posVelToOrbit(cls, x, v, t = 0.0, att_ = defaultAttractor):
        orbit_ = cls()
        _, orbitDict, _ = oc.posVelToStruct(x, v,  att_.mu, t, nout = 3)
        orbit_.p    = orbitDict['p']
        orbit_.e    = orbitDict['e']
        orbit_.th0  = orbitDict['aop']
        orbit_.incl = orbitDict['incl']
        orbit_.lan  = orbitDict['lan']
        orbit_.t0   = orbitDict['peT']
        orbit_.att  = att_
        return orbit_
    
    @classmethod
    def sfsDictToOrbit(cls, orbitDictionary, attractorDictionary = None, attractor = defaultAttractor):
        # Get eccentricity with parabolic orbits corrected as hyperbolic
        e_ = float(orbitDictionary['ECC'])
        if e_ == 1.0:
            e_ += np.exp2(-52)
        # Get orbit plane angles
        th0_  = float(orbitDictionary['LPE'])*degToRad
        incl_ = float(orbitDictionary['INC'])*degToRad
        lan_  = float(orbitDictionary['LAN'])*degToRad
        # Additional quantities to be processed
        semiMajorAxis = float(orbitDictionary['SMA'])
        referenceTime = float(orbitDictionary['EPH'])
        meanAnomaly   = float(orbitDictionary['MNA'])
        attractorID   = int(orbitDictionary['REF'])
        # Get attractor
        if len(cls.system) > 0 and attractorID in cls.system:
            att_ = cls.system[attractorID]
        elif attractorDictionary != None and attractorID in attractorDictionary:
            att_ = attractorDictionary[attractorID]
        else:
            att_ = attractor
        # Get semi-latus rectum
        p_ = semiMajorAxis*(1.0 - e_*e_)
        # Get reference time at periapsis
        t0_ = referenceTime - (meanAnomaly*np.sqrt((semiMajorAxis**3)/att_.mu))
        return cls.defineParameters(p_, e_, th0_, incl_, lan_, t0_, att_)
    
    @classmethod
    def jsonDictToOrbit(cls, objectDictionary, attractor = defaultAttractor):
        # Get eccentricity with parabolic orbits corrected as hyperbolic
        e_ = objectDictionary['orbitalCharacteristics']['orbitalEccentricity']
        if e_ == 1.0:
            e_ += np.exp2(-52)
        # Get orbit plane angles
        th0_  = objectDictionary['orbitalCharacteristics']['argumentOfPeriapsis']*degToRad
        incl_ = objectDictionary['orbitalCharacteristics']['orbitalInclination']*degToRad
        lan_  = objectDictionary['orbitalCharacteristics']['longitudeOfTheAscendingNode']*degToRad
        # Additional quantities to be processed
        semiMajorAxis = objectDictionary['orbitalCharacteristics']['semimajorAxis']
        referenceTime = 0.0
        meanAnomaly   = objectDictionary['orbitalCharacteristics']['meanAnomaly']
        attractorID   = objectDictionary['attractor']['REF']
        # Get attractor
        if len(cls.system) > 0 and attractorID in cls.system:
            att_ = cls.system[attractorID]
        else:
            att_ = attractor
        # Get semi-latus rectum
        p_ = semiMajorAxis*(1.0 - e_*e_)
        # Get reference time at periapsis
        t0_ = referenceTime - (meanAnomaly*np.sqrt((semiMajorAxis**3)/att_.mu))
        return cls.defineParameters(p_, e_, th0_, incl_, lan_, t0_, att_)

    def orbitToSfsDict(self, attractorDictionary = None, attractorID_ = 1):
        # Get attractor
        if len(self.system) > 0 and self.att in list(self.system.values()):
            attractorID = list(self.system.keys())[list(self.system.values()).index(self.att)]
        elif attractorDictionary != None and self.att in list(attractorDictionary.values()):
            attractorID = list(attractorDictionary.keys())[list(attractorDictionary.values()).index(self.att)]
        else:
            attractorID = attractorID_
        # Additional quantities to be processed
        semiMajorAxis = self.semiMajorAxis()
        referenceTime = self.t0
        meanAnomaly   = 0.0
        # Get orbit plane angles
        argumentOfPeriapsis       = radToDeg*self.th0
        orbitalInclination        = radToDeg*self.incl
        longitudeOfAscendingNode  = radToDeg*self.lan
        # Initialize dictionary    
        orbitDictionary = {}
        # Get eccentricity with parabolic orbits corrected as hyperbolic    
        eccentricity = self.e
        if eccentricity == 1.0:
            eccentricity += np.exp2(-52)
        orbitDictionary['SMA'] = f"{semiMajorAxis:.17g}"
        orbitDictionary['ECC'] = f"{eccentricity:.17g}"
        orbitDictionary['INC'] = f"{orbitalInclination:.17g}"
        orbitDictionary['LPE'] = f"{argumentOfPeriapsis:.17g}"
        orbitDictionary['LAN'] = f"{longitudeOfAscendingNode:.17g}"
        orbitDictionary['MNA'] = f"{meanAnomaly:.17g}"
        orbitDictionary['EPH'] = f"{referenceTime:.17g}"
        orbitDictionary['REF'] = f"{attractorID:d}"
        return orbitDictionary
    
    def matlabDict(self):
        orbit_ = { 'p'    : self.p, \
                   'e'    : self.e, \
                   'lan'  : self.lan, \
                   'incl' : self.incl, \
                   'aop'  : self.th0, \
                   'peT'  : self.t0}
        return orbit_

    def timeToPosVel(self, t):
        x, v, _ = oc.structTimeToPosVelArg(t, self.matlabDict(), self.att.mu, nout = 3)
        return x, v
    
    def arumentToPosVel(self, th):
        x, v, _ = oc.structArgToPosVelTime(th, self.matlabDict(), self.att.mu, nout = 3)
        return x, v
    
    def refFrameMatrix(self):
        return oc.rOrb(self.lan, self.incl, self.th0)
    
    def argumentToTime(self, th):
        return oc.argToTime(th, self.p, self.e, self.att.mu, self.t0, nout = 1)
    
    def timeToArgument(self, t):
        return oc.timeToArg(t, self.p, self.e, self.att.mu, self.t0, nout = 1)
    
    def angularVelocity(self, th):
        return oc.angularVelocity(th, self.p, self.e, self.att.mu, nout = 1)

    def periapsis(self):
        return self.p / (1.0 + self.e)
    
    def apoapsis(self):
        return self.p / max(1.0 - self.e, tiny*self.p)

    def semiMajorAxis(self):
        return self.p / max(1.0 - self.e**2, tiny*self.p)

    def semiMinorAxis(self):
        return self.p / max(np.sqrt(1.0 - self.e**2), tiny*self.p)