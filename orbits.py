import numpy as np
import oct2py as op
from oct2py import octave as oc
"""
Provides an interface between python and the orbital mechanics library written 
by me for matlab a few years ago

Created on Tue Nov 12 15:18:07 2024

@author: Tommaso Zanelli
"""

path=r"OctaveScripts\\"

oc.addpath(path)
oc.addpath(oc.genpath(path))

tiny = 5.56268464626801e-309

class attractor:
    def __init__(self, mu_ = 1.0, name_=""):
        self.mu   = mu_
        self.name = name_
        
defaultAtt = attractor()
        
class orbit:
    def __init__(self):
        self.p    = 1.0
        self.e    = 0.0
        self.th0  = 0.0
        self.incl = 0.0
        self.lan  = 0.0
        self.t0   = 0.0
        self.att  = defaultAtt

    @classmethod
    def defineParameters(cls, att_, p_, e_ = 0.0, th0_ = 0.0, incl_ = 0.0, lan_ = 0.0, t0_ = 0.0):
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
    def posVelToOrbit(cls, att_, x, v, t=0.0):
        orbit_ = cls()
        _, orbitDict, _ = oc.posVelToStruct(x, v,  att_.mu, t, nout=3)
        orbit_.p    = orbitDict['p']
        orbit_.e    = orbitDict['e']
        orbit_.th0  = orbitDict['aop']
        orbit_.incl = orbitDict['incl']
        orbit_.lan  = orbitDict['lan']
        orbit_.t0   = orbitDict['peT']
        orbit_.att  = att_
        return orbit_
    
    def matlabDict(self):
        orbit_ = { 'p'    : self.p, \
                   'e'    : self.e, \
                   'lan'  : self.lan, \
                   'incl' : self.incl, \
                   'aop'  : self.th0, \
                   'peT'  : self.t0}
        return orbit_

    def timeToPosVel(self, t):
        x, v, _ = oc.structTimeToPosVelArg(t, self.matlabDict(), self.att.mu, nout=3)
        return x, v
        
    def arumentToPosVel(self, th):
        x, v, _ = oc.structArgToPosVelTime(th, self.matlabDict(), self.att.mu, nout=3)
        return x, v
    
    def refFrameMatrix(self):
        return oc.rOrb(self.lan, self.incl, self.th0)
    
    def argumentToTime(self, th):
        return oc.argToTime(th, self.p, self.e, self.att.mu, self.t0, nout=1)
        
    def timeToArgument(self, t):
        return oc.timeToArg(t, self.p, self.e, self.att.mu, self.t0, nout=1)
    
    def angularVelocity(self, th):
        return oc.angularVelocity(th, self.p, self.e, self.att.mu, nout=1)

    def periApsis(self):
        return self.p / (1.0 + self.e)
    
    def apoApsis(self):
        return self.p / max(1.0 - self.e, tiny*self.p)

    def sMajAxis(self):
        return self.p / max(1.0 - self.e**2, tiny*self.p)

    def sMinAxis(self):
        return self.p / max(np.sqrt(1.0 - self.e**2), tiny*self.p)

class celestialBody(attractor):
    def __init__(self, mu_ = 1.0, name_="", orbit_=None, rSOI_=None, satList_ = []):
        self.mu      = mu_
        self.name    = name_
        self.orbit   = orbit_
        self.rSOI    = rSOI_
        self.satList = satList_
        if self.rSOI == None:
            self.rSOI = 1.0/tiny
            if self.orbit != None:
                self.rSOI = self.computeSOI()
    
    def computeSOI(self):
        mG = self.mu
        MG = self.orbit.att.mu
        a  = self.orbit.sMajAxis()
        return a*np.pow2(mG/MG, 0.4)

    def addSatellite(self, sat):
        self.satList.append(sat)

