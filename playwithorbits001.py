#!/usr/bin/env python3
import pyksp as pk
import numpy as np
import os
import pickle
import tqdm

from optiKev import optimizeCost

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.close('all')

tck = pk.timer()
tck2 = pk.timer()
tck.tic()
saveDict = pk.sfsFileToDict(r"C:\GOG\Kerbal Space Program\saves\Career_01\persistent.sfs")

time = float(saveDict['GAME']['FLIGHTSTATE']['UT'])
dt = 3600*6*40 # 30 Kerbin days

oK = pk.orbit.system[1].orbit

kevinId = [vessel['name'] for vessel in saveDict['GAME']['FLIGHTSTATE']['VESSEL']].index('Kevin I')
lrId = [vessel['name'] for vessel in saveDict['GAME']['FLIGHTSTATE']['VESSEL']].index('Lexa Raider')

oKev = pk.orbit.sfsDictToOrbit(saveDict['GAME']['FLIGHTSTATE']['VESSEL'][kevinId]['ORBIT'])
oLr =  pk.orbit.sfsDictToOrbit(saveDict['GAME']['FLIGHTSTATE']['VESSEL'][lrId]['ORBIT'])

# When does the asteroid enter Kerbin orbit?
tEntry, d, th1, th2 = pk.entryInSphere(time,  time + dt, oLr, pk.orbit.system[1])


xLrEntry, vLrEntry = oLr.timeToPosVel(tEntry)
xKEntry, vKEntry = oK.timeToPosVel(tEntry)

# Orbit after Kerbin entry
oLrK = pk.orbit.posVelToOrbit(xLrEntry - xKEntry, vLrEntry - vKEntry, t = tEntry, att_ = pk.orbit.system[1])

# Kerbin exit time
tExit = oLrK.argumentToTime(-oLrK.timeToArgument(tEntry))

xLrExit, vLrExit = oLrK.timeToPosVel(tExit)
xKExit, vKExit = oK.timeToPosVel(tExit)

# Orbiit after exit
oLr2 = pk.orbit.posVelToOrbit(xLrExit + xKExit, vLrExit + vKExit, t = tExit, att_ = pk.orbit.system[0])

tt1 = np.linspace(time, tEntry, 5001)
tt2 = np.linspace(tEntry, tExit, 5001)
tt3 = np.linspace(tExit, time + dt, 5001)

xKerbin1, vKerbin1 = oK.timeToPosVel(tt1)
xKerbin2, vKerbin2 = oK.timeToPosVel(tt2)
xKerbin3, vKerbin3 = oK.timeToPosVel(tt3)

xKerbin = np.concatenate((xKerbin1, xKerbin2, xKerbin3), 1)
vKerbin = np.concatenate((vKerbin1, vKerbin2, vKerbin3), 1)

xLr1, _ = oLr.timeToPosVel(tt1)
xLr2, _ = oLrK.timeToPosVel(tt2)
xLr3, _ = oLr2.timeToPosVel(tt3)
tt = np.concatenate((tt1, tt2, tt3))

xLr2_ = xLr2 + xKerbin2

_xKev, _ = oKev.timeToPosVel(tt)
xKev = _xKev + xKerbin


def timeToPosVel_Lr(t):
    scalar = np.ndim(t) == 0
    t = np.atleast_1d(np.asarray(t, dtype=float))
    
    x = np.zeros((3, len(t)))
    v = np.zeros((3, len(t)))
    
    m1 = t < tEntry
    m2 = (t >= tEntry) & (t <= tExit)
    m3 = t > tExit
    
    if m1.any():
        x[:, m1], v[:, m1] = oLr.timeToPosVel(t[m1])
    if m2.any():
        xKerbin, vKerbin = oK.timeToPosVel(t[m2])
        x[:, m2], v[:, m2] = oLrK.timeToPosVel(t[m2])
        x[:, m2] = x[:, m2] + xKerbin
        v[:, m2] = v[:, m2] + vKerbin
    if m3.any():
        x[:, m3], v[:, m3] = oLr2.timeToPosVel(t[m3])
    
    if scalar:
        return x[:, 0], v[:, 0]
    return x, v

def ff(t):
    return oLrK.timeToPosVel(t)[0][2]

# Time at which LR is co-planar with Kerbin
tMax = 5355865.8685875

t1 = time + 0.2*(tMax - time)
t2 = time + 0.8*(tMax - time)

dv = np.array([-939.50512533,   854.34996514, -2775.7143412])

def timeToPosVel_Kev(t, t1, dv):
    scalar = np.ndim(t) == 0
    t = np.atleast_1d(np.asarray(t, dtype=float))
    
    x = np.zeros((3, len(t)))
    v = np.zeros((3, len(t)))

    oKevM1 = pk.applyManoeuvre(oKev, t1, dv)
    if oKevM1.e > 1.0 or oKevM1.apoapsis() > oKevM1.att.rSOI:
        tExitKev = oKevM1.argumentToTime(np.acos(((oKevM1.p/oKevM1.att.rSOI) - 1)/(oKevM1.e)))
        xExitKev, vExitKev = oKevM1.timeToPosVel(tExitKev)
        xExitKevKer, vExitKevKer = oK.timeToPosVel(tExitKev)
        oKevM2 = pk.orbit.posVelToOrbit(xExitKev + xExitKevKer, vExitKev + vExitKevKer, t = tExitKev, att_ = oK.att)
    else:
        tExitKev = np.exp2(1023)
        oKevM2 = oKevM1

    m1 = t < t1
    m2 = (t >= t1) & (t <= tExitKev)
    m3 = t > tExitKev
    
    if m1.any():
        xKerbin1, vKerbin1 = oK.timeToPosVel(t[m1])
        x[:, m1], v[:, m1] = oKev.timeToPosVel(t[m1])
        x[:, m1] = x[:, m1] + xKerbin1
        v[:, m1] = v[:, m1] + vKerbin1
    if m2.any():
        xKerbin2, vKerbin2 = oK.timeToPosVel(t[m2])
        x[:, m2], v[:, m2] = oKevM1.timeToPosVel(t[m2])
        x[:, m2] = x[:, m2] + xKerbin2
        v[:, m2] = v[:, m2] + vKerbin2
    if m3.any():
        x[:, m3], v[:, m3] = oKevM2.timeToPosVel(t[m3])
    
    if scalar:
        return x[:, 0], v[:, 0], tExitKev, oKevM1, oKevM2
    return x, v, tExitKev, oKevM1, oKevM2

_, _, tExitKev, oKevM1, oKevM2 = timeToPosVel_Kev(t2, t1, dv)

tKm = np.linspace(t1, time + dt, 15001)
tKm1 = tKm[np.where(tKm <= tExitKev)[0]] 
tKm2 = tKm[np.where(tKm > tExitKev)[0]] 

_xKevPM1, _ = oKevM1.timeToPosVel(tKm1)
xKerbinKev, _ = oK.timeToPosVel(tKm1)
xKevPM1 = _xKevPM1 + xKerbinKev
xKevPM2, _ = oKevM2.timeToPosVel(tKm2)

xKevMM, _, tExitKev_ck, oKevM1_ck, oKevM2_ck = timeToPosVel_Kev(tt, t1, dv)

def distVel(t1, t2, dv):
    xKev, vKev, _, _, _ = timeToPosVel_Kev(t2, t1, dv)
    xLr, vLr = timeToPosVel_Lr(t2)
    return np.linalg.norm(xKev - xLr), np.linalg.norm(vKev - vLr)


xLr2ck, vLr2ck = timeToPosVel_Lr(tt)

def xvJacobian(dv, t1, t2, ddv = []):
    if len(ddv) == 0:
        ddv = np.array([0,0,0])
    ddv = np.maximum(ddv, np.exp2(np.ceil(np.log2(np.maximum(abs(dv), 1)) - 26)))
    J = np.zeros((3, 3))
    x, _, _, _, _ = timeToPosVel_Kev(t2, t1, dv)
    for ii in range(3):
        ddv2 = np.zeros(3)
        ddv2[ii] = ddv[ii]
        xp, _, _, _, _ = timeToPosVel_Kev(t2, t1, dv + ddv2)
        xm, _, _, _, _ = timeToPosVel_Kev(t2, t1, dv - ddv2)
        for jj in range(3):
            J[ii][jj] = 0.5*(xp[jj] - xm[jj])/ddv[ii]
    return J

# def optimizeGivenTime(t1, t2, dv0 = np.array([0.0, 0.0, 0.0]), verbose = False, errMax = None):
#     dv2 = dv0.copy()
#     iterate = True
#     xLr, _ = timeToPosVel_Lr(t2)
#     if errMax == None:
#         errMax = np.exp2(np.log2(np.linalg.norm(xLr)) - 45)
#     while iterate:
#         x, _, _, _, _ = timeToPosVel_Kev(t2, t1, dv2)
#         dx = xLr - x
#         C0 = np.linalg.norm(dx)
#         if (verbose): print(C0)
#         J = xvJacobian(dv2, t1, t2)
#         dva = np.linalg.solve(J.transpose(), dx)
#         a = 1
#         while distVel(t1, t2, dv2 + a*dva)[0] > C0 and a > 1e-8 :
#          	a /= 2.0
#         dv2 = dv2 + a*dva
#         C1 = distVel(t1, t2, dv2)[0]
#         if (verbose): print(C1, a, '\n')
#         if C1 >= C0 or C1 <= errMax:
#             iterate = False
#     return dv2

def optimizeGivenTime(t1, t2, dv0 = np.array([0.0, 0.0, 0.0]), verbose = False, errMax = None, xLr = [], maxIter = None):
    dv2 = dv0.copy()
    iterate = True
    if len(xLr) == 0: 
        xLr, _ = timeToPosVel_Lr(t2)
    if errMax == None:
        errMax = np.exp2(np.log2(np.linalg.norm(xLr)) - 45)
    if maxIter == None:
        maxIter = 2147483648
    nIter = 0
    while iterate:
        x, _, _, _, _ = timeToPosVel_Kev(t2, t1, dv2)
        dx = xLr - x
        C0 = np.linalg.norm(dx)
        #if (verbose):
        #    print("                                          ", end = '\r')
        #    print(C0, end='\r')
        ddv = np.exp2(np.ceil(np.log2(np.maximum(abs(dv2), 1)) - 26))
        J = np.zeros((3, 3))
        HJ = np.zeros((3, 3, 3))
        ddv3 = np.zeros(3)
        for ii in range(3):
            ddv3[ii] = ddv[ii]
            xp, _, _, _, _ = timeToPosVel_Kev(t2, t1, dv2 + ddv3)
            xm, _, _, _, _ = timeToPosVel_Kev(t2, t1, dv2 - ddv3)
            for jj in range(3):
                J[jj][ii] = 0.5*(xp[jj] - xm[jj])/ddv[ii]
            ddv3[ii] = 0.0
        dva = np.linalg.solve(J, dx)
        dva2 = dva.copy()
        a = 1.0
        xKev, vKev, _, _, _ = timeToPosVel_Kev(t2, t1, dv2 + dva2)
        C1 = np.linalg.norm(xKev - xLr)
        while C1 > C0 and a > 1e-8 and not np.any(np.isnan(dva2)):
            a /= 2.0
            dva2 /= 2.0
            xKev, vKev, _, _, _ = timeToPosVel_Kev(t2, t1, dv2 + dva2)
            C1 = np.linalg.norm(xKev - xLr)
        dv2 = dv2 + dva2
        if (verbose): 
            print("                                          ", end = '\r')
            print(nIter, C1, a, end = '\r')
        nIter += 1
        if C1 >= C0 or C1 <= errMax or nIter > maxIter or np.any(np.isnan(dv2)):
            iterate = False
    return dv2, xKev, vKev


dv2, _, _ = optimizeGivenTime(t1, t2, dv, verbose = True)

xLr, vLr = timeToPosVel_Lr(t2)
xKevO, vKevO, _, _, _ = timeToPosVel_Kev(t2, t1, dv2)

print("Final distance: ", np.linalg.norm(xLr - xKevO), " [m]")
dv1n = np.linalg.norm(dv2)
dv2n = np.linalg.norm(vLr - vKevO)
print("dv1 = ", dv1n, " [m/s];  dv2 = ", dv2n, " [m/s], tot : ", dv1n + dv2n, " [m/s]")

def costTime(t1, t2, dv0 = np.array([0.0, 0.0, 0.0]), verbose = False, errMax = None, maxIter = None):
    xLr, vLr = timeToPosVel_Lr(t2)
    dv2, xKevO, vKevO = optimizeGivenTime(t1, t2, dv0, verbose = verbose, errMax = errMax, xLr = xLr, maxIter = maxIter)
    dv1n = np.linalg.norm(dv2)
    dv2n = np.linalg.norm(vLr - vKevO)
    return dv1n + dv2n, dv2, np.linalg.norm(xLr - xKevO)

def optimiset2(t1, t20 = None, dv0 = np.array([0.0, 0.0, 0.0]), verbose = False, maxDv = None, innerErrMax = 100, outerErrMax = None):
    t2l = [t1 + 2000, tMax - 2000]
    iterate = True
    if t20 == None:
        t2 = (t2l[0] + t2l[1])*0.5
    else:
        t2 = t20 
    if verbose: print("Computing initial cost:")
    C0, dv, _ = costTime(t1, t2, dv0 = dv0, verbose = verbose, errMax = innerErrMax)
    if maxDv == None:
        maxDv = np.exp2(np.log2(C0) - 26)
    if verbose: print("Initial cost is: ", C0)
    while iterate:
        ddt = np.exp2(np.ceil(np.log2(np.maximum(abs(dt), 1)) - 26))
        t2p =  min(t2 + ddt, t2l[1])
        if verbose: print("Computing with increased t2: ", t2p)
        Cp, _, _ = costTime(t1, t2p, dv0 = dv, verbose = verbose, errMax = innerErrMax)
        if verbose: print("Increased t2 cost is: ", Cp)
        
        t2m = max(t2 - ddt, t2l[0])
        if verbose: print("Computing with decreased t2: ", t2m)
        Cm, _, _ = costTime(t1, t2m, dv0 = dv, verbose = verbose, errMax = innerErrMax)
        if verbose: print("Decreased t2 cost is: ", Cm)
        dC = 0.5*(Cp - Cm)/ddt
        dC2 = (Cp + Cm - 2.0*C0)/(ddt**2)
        dt2 = -dC/dC2
        if verbose: print("Newton dt: ", dt2)
        a = 1
        if verbose: print("Iterating while a = ", a)
        CC, dvn, _ = costTime(t1, min(max(t2 + a*dt2, t2l[0]), t2l[1]), dv0 = dv, verbose = verbose, errMax = innerErrMax)
        if verbose: print("Cost is = ", CC)
        while CC > C0 and a > 1.220703125e-4:
            a /= 2.0
            if verbose: print("CC = ", CC, "Iterating while a = ", a)
            CC, dvn, _ = costTime(t1, min(max(t2 + a*dt2, t2l[0]), t2l[1]), dv0 = dvn, verbose = verbose, errMax = innerErrMax)   
            if verbose: print("Cost is = ", CC)
        t2 = min(max(t2 + a*dt2, t2l[0]), t2l[1])
        if verbose: print("Computing updated cost:")
        C1, dv, _ = costTime(t1, t2, dv0 = dvn, verbose = verbose, errMax = outerErrMax)
        if C1 >= C0 or abs(C1 - C0) < maxDv:
            iterate = False
        if (verbose): print("Best dv so far: ", C1, a, '\n')
        C0 = C1
    return t2


# tck2.tic()
# t2n = optimiset2(t1, t2, dv2, verbose = True, maxDv = 1, innerErrMax = 500, outerErrMax = 0.1)
# tck2.toc("Optimization time: ")
# print("Previous optimal t2 = ", t2, "[s], Previous optimum: ", costTime(t1, t2, dv0 = dv2, verbose = False), " [m/s]")
# print("Current  optimal t2 = ", t2n, "[s], Current  optimum: ", costTime(t1, t2n, dv0 = dv2, verbose = False), " [m/s]")


TKev = 2*np.pi*np.sqrt((oKev.semiMajorAxis()**3)/oKev.att.mu)
t1Lims = [time + 2000, tMax - 4000]
t1v = np.linspace(t1Lims[0], t1Lims[1], 8*int(np.ceil(np.diff(t1Lims)[0]/(TKev))) - 3)
t2v = t1v.copy() + 2000

#dvt1v0 = np.array([-1083.77897426, 2824.18757739, -860.81635755])
dvt1v0 = np.array([-5405.33668948,   172.65934297,  -369.3795052 ])

dt1 = np.average(np.diff(t1v))

CostFile = "costFile.pkl"
print('\n\n')
if os.path.exists(CostFile):
    fileId = open(CostFile, 'rb')
    dvl = pickle.load(fileId)
    fileId.close()
else:
    vvv = False
    #dvn = dvt1v0.copy()
    t2n = t2
    dvl = []
    #dvll = []
    dtv = np.linspace(0, dt1, 64)
    dtdt = np.average(np.diff(dtv))
    for ii, t1tc in enumerate(t1v): 
        t2tc = t2v[-1]
        tck2.tic()
        C = 0
        if len(dvl) == 0:
            dvn0 = dvt1v0.copy() 
        elif len(dvl) == 1:
            dvn0 = np.array([0.0, 0.0, 0.0])
        elif len(dvl) == 2:
            dvn0[0] = np.polyval(np.polyfit([t[0][0] for t in dvl[-2:]], [vv[0][2][0] for vv in dvl[-2:]], 1), t1tc)
            dvn0[1] = np.polyval(np.polyfit([t[0][0] for t in dvl[-2:]], [vv[0][2][1] for vv in dvl[-2:]], 1), t1tc)
            dvn0[2] = np.polyval(np.polyfit([t[0][0] for t in dvl[-2:]], [vv[0][2][2] for vv in dvl[-2:]], 1), t1tc)
        else:
            dvn0[0] = np.polyval(np.polyfit([t[0][0] for t in dvl[-3:]], [vv[0][2][0] for vv in dvl[-3:]], 2), t1tc)
            dvn0[1] = np.polyval(np.polyfit([t[0][0] for t in dvl[-3:]], [vv[0][2][1] for vv in dvl[-3:]], 2), t1tc)
            dvn0[2] = np.polyval(np.polyfit([t[0][0] for t in dvl[-3:]], [vv[0][2][2] for vv in dvl[-3:]], 2), t1tc)
        C, dvn, D = costTime(t1tc, t2tc, dv0 = dvn0, verbose = True, maxIter = 128)
        if D < 100 and not np.isnan(D):
            dvl.append([[t1tc, t2tc, dvn, C, D]])
            print("t1 = ", t1tc, "[s]; t2 = ", t2tc, "[s]; dv = ", C, "[m/s], dist. = ", D, " [m]", end = '\n')
        else:
            print("CONVERGENCE FAILURE: D = ", D, " [m]", end = '\n')
        tck2.toc("  Optimization time: ")
    fileId = open(CostFile, 'wb')
    pickle.dump(dvl, fileId)
    fileId.close()

addendum_1 = False
if addendum_1:
    dvl2 = dvl.copy()
    vjj = [np.zeros(3) for i in range(4)]
    dvn = np.zeros(3)
    t2tc = t2v[-1]
    xt, _ = timeToPosVel_Lr(t2tc)
    for ii, t1tc in enumerate(t1v):
        tts = np.array([t[0][0] for t in dvl])
        diffts = tts - t1tc
        if min(abs(diffts)) < 0.01:
            print(f"{t1tc:g} : FOUND!")
        else:
            jj = np.where(np.diff(np.sign(diffts)) != 0)[0][0]
            dvnn = dvl[jj][0][2].copy()
            if jj == 0:
                jv = [jj, jj + 1, jj + 2, jj + 3]
            elif jj == len(dvl) - 2:
                jv = [jj - 2, jj - 1, jj, jj + 1]
            elif jj == len(dvl) - 1:
                jv = [jj - 3, jj - 2, jj - 1, jj]
            else:
                jv = [jj - 1, jj, jj + 1, jj + 2]
            print(f"{t1tc:g} Not found. Computing...")
            for j, jj in enumerate(jv):
                xj, vj, _, _, _ = timeToPosVel_Kev(tts[jj], tMax, np.zeros(3))
                Tj = pk.localReferenceFrame(np.reshape(xj, (3, 1)), np.reshape(vj, (3, 1)))
                vjj[j] = Tj@dvl[jj][0][2]  + vj - (xt - xj)/(t2tc - tts[jj])
            xm, vm, _, _, _ = timeToPosVel_Kev(t1tc, tMax, np.zeros(3))
            Tm = pk.localReferenceFrame(np.reshape(xm, (3, 1)), np.reshape(vm, (3, 1)))
            for k in range(3):
                dvn[k] = np.polyval(np.polyfit(tts[jv], [q[k] for q in vjj], 3), t1tc)
            dvn = Tm.transpose()@(dvn - vm + (xt - xm)/(t2tc - t1tc))
            C, dvnn, D = costTime(t1tc, t2tc, dv0 = dvn, verbose = True, maxIter = 48)
            if D >= 100 or np.isnan(D):
                C, dvnn, D = costTime(t1tc, t2tc, dv0 = dvnn, verbose = True, maxIter = 48)
            if D >= 100 or np.isnan(D):
                C, dvnn, D = costTime(t1tc, t2tc, dv0 = 0.5*(dvl[jj][0][2] + dvl[jj + 1][0][2]), verbose = True, maxIter = 48)
            if D >= 100 or np.isnan(D):
                C, dvnn, D = costTime(t1tc, t2tc, dv0 = np.zeros(3), verbose = True, maxIter = 48)
            if D < 100 and not np.isnan(D):
                dvl2.append([[t1tc, t2tc, dvnn, C, D]])
                print("t1 = ", t1tc, "[s]; t2 = ", t2tc, "[s]; dv = ", C, "[m/s], dist. = ", D, " [m]", end = '\n')
            else:
                print("CONVERGENCE FAILURE: D = ", D, " [m]", end = '\n')
    dvl2.sort(key = lambda x: x[0][0])
    dvl = dvl2
    fileId = open(CostFile, 'wb')
    pickle.dump(dvl, fileId)
    fileId.close()


dvl2 = dvl.copy()
vjj = [np.zeros(3) for i in range(4)]
dvn = np.zeros(3)
t2tc = t2v[-1]
xt, _ = timeToPosVel_Lr(t2tc)
for ii, t1tc in enumerate(t1v[8:9]): # 11
    tts = np.array([t[0][0] for t in dvl])
    diffts = tts - t1tc
    if min(abs(diffts)) < 0.01:
        jj = np.where(abs(diffts) == min(abs(diffts)))[0][0]
        print(f"{t1tc:g} : FOUND!")
        t1tc = dvl[jj][0][0] 
        C = dvl[jj][0][3] 
        D = dvl[jj][0][4] 
        print("t1 = ", t1tc, "[s]; t2 = ", t2tc, "[s]; dv = ", C, "[m/s], dist. = ", D, " [m]", end = '\n')
        t1n, t2n, dvnn, Ci, Cf = optimizeCost(costTime, t1tc, t2tc, dvl[jj][0][2], verbose = True, maxIter = 48)
        print("t1 = ", t1n, "[s]; t2 = ", t2n, "[s]; dv = ", Cf, "[m/s], dist. = ", 1, " [m]", end = '\n')
    else:
        print(f"{t1tc:g} Not found. Computing...")
        # jj = np.where(np.diff(np.sign(diffts)) != 0)[0][0]
        # dv01 = dvl[jj][0][2]
        # t1 = dvl[jj][0][0]
        # dt1 = np.exp2(np.floor(np.log2(abs(t1))) - 21)
        # dt1m1 = 1.0/dt1
        # C01, dv01, D01 = costTime(t1, t2tc, dv0 = dv01, verbose = True)
        # Cp1, dvp1, Dp1 = costTime(t1 + dt1, t2tc, dv0 = dv01, verbose = True)
        # Cm1, dvm1, Dm1 = costTime(t1 - dt1, t2tc, dv0 = dv01, verbose = True)
        # dvt1 = dv01 + 0.5*(dvp1 - dvm1)*(t1tc - t1)*dt1m1 + 0.5*(dvp1 - 2.0*dv01 + dvm1)*(t1tc - t1)*(t1tc - t1)*dt1m1*dt1m1
        # dv02 = dvl[jj + 1][0][2]
        # t2 = dvl[jj + 1][0][0]
        # dt2 = np.exp2(np.floor(np.log2(abs(t2))) - 21)
        # dt2m1 = 1.0/dt2
        # C02, dv02, D02 = costTime(t2, t2tc, dv0 = dv02, verbose = True)
        # Cp2, dvp2, Dp2 = costTime(t2 + dt2, t2tc, dv0 = dv02, verbose = True)
        # Cm2, dvm2, Dm2 = costTime(t2 - dt2, t2tc, dv0 = dv02, verbose = True)
        # dvt2 = dv02 + 0.5*(dvp2 - dvm2)*(t2tc - t2)*dt2m1 + 0.5*(dvp2 - 2.0*dv02 + dvm2)*(t1tc - t2)*(t1tc - t2)*dt2m1*dt2m1

        # for j, jj in enumerate(jv):
        #     xj, vj, _, _, _ = timeToPosVel_Kev(tts[jj], tMax, np.zeros(3))
        #     Tj = pk.localReferenceFrame(np.reshape(xj, (3, 1)), np.reshape(vj, (3, 1)))
        #     vjj[j] = Tj@dvl[jj][0][2]  + vj - (xt - xj)/(t2tc - tts[jj])
        # xm, vm, _, _, _ = timeToPosVel_Kev(t1tc, tMax, np.zeros(3))
        # Tm = pk.localReferenceFrame(np.reshape(xm, (3, 1)), np.reshape(vm, (3, 1)))
        # for k in range(3):
        #     dvn[k] = np.polyval(np.polyfit(tts[jv], [q[k] for q in vjj], 3), t1tc)
        # dvn = Tm.transpose()@(dvn - vm + (xt - xm)/(t2tc - t1tc))
        # C, dvnn, D = costTime(t1tc, t2tc, dv0 = dvn, verbose = True, maxIter = 48)
        # if D >= 100 or np.isnan(D):
        #     C, dvnn, D = costTime(t1tc, t2tc, dv0 = dvnn, verbose = True, maxIter = 48)
        # if D >= 100 or np.isnan(D):
        #     C, dvnn, D = costTime(t1tc, t2tc, dv0 = 0.5*(dvl[jj][0][2] + dvl[jj + 1][0][2]), verbose = True, maxIter = 48)
        # if D >= 100 or np.isnan(D):
        #     C, dvnn, D = costTime(t1tc, t2tc, dv0 = np.zeros(3), verbose = True, maxIter = 48)
        # if D < 100 and not np.isnan(D):
        #     dvl2.append([[t1tc, t2tc, dvnn, C, D]])
        #     print("t1 = ", t1tc, "[s]; t2 = ", t2tc, "[s]; dv = ", C, "[m/s], dist. = ", D, " [m]", end = '\n')
        # else:
        #     print("CONVERGENCE FAILURE: D = ", D, " [m]", end = '\n')        t1tc = dvl[jj][0][0] 
        C = None #dvl[jj][0][3] 
        D = None #dvl[jj][0][4] 
        print("t1 = ", t1tc, "[s]; t2 = ", t2tc, "[s]; dv = ", C, "[m/s], dist. = ", D, " [m]", end = '\n')
dvl2.sort(key = lambda x: x[0][0])
dvl = dvl2
#fileId = open(CostFile, 'wb')
#pickle.dump(dvl, fileId)
#fileId.close()

# CostFile = "costFile.pkl"
# if os.path.exists(CostFile):
#     fileId = open(CostFile, 'rb')
#     dvl = pickle.load(fileId)
#     fileId.close()
# else:
#     print('\n\n')
#     dvn = dv2.copy()
#     t2n = t2
#     dvl = []
#     dvll = []
#     for t1tc in t1v:
#         if len(dvl) == 0:
#             dvn = dvt1v0.copy()
#         elif len(dvl) == 1:
#             dvn = dvl[0][0][2]
#         else:
#             dvn = 2.0*dvl[-1][0][2] - dvl[-2][0][2]
#         dvll.clear()
#         for t2tc in np.flip(t2v[np.where(t2v > t1tc + 2000)[0]]):
#             if len(dvll) > 1:
#                 dvn = 2.0*dvll[-1][2] - dvll[-2][2]
#             tck2.tic()
#             C = 0
#             C, dvn, D = costTime(t1tc, t2tc, dv0 = dvn, verbose = False)
#             dvll.append([t1tc, t2tc, dvn, C, D])
#             print("t1 = ", t1tc, "[s]; t2 = ", t2tc, "[s]; dv = ", C, "[m/s], dist. = ", D, " [m]", end = '')
#             tck2.toc("  Optimization time: ")
#         dvl.append(dvll.copy())
#     fileId = open(CostFile, 'wb')
#     pickle.dump(dvl, fileId)
#     fileId.close()

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

ax.plot(xLr1[0], xLr1[1], xLr1[2], color = [0, 0, 1])
ax.plot(xLr2_[0], xLr2_[1], xLr2_[2], color = [1, 0, 0])
ax.plot(xLr3[0], xLr3[1], xLr3[2], color = [0, 0, 1])
ax.plot(xKev[0], xKev[1], xKev[2], color = [0, 1, 0])
ax.plot(xLr2ck[0], xLr2ck[1], xLr2ck[2], '--', color = [0, 0, 0])
ax.plot(xKevPM1[0], xKevPM1[1], xKevPM1[2], color = [0.5, 0.5, 0], linewidth = 3)
ax.plot(xKevPM2[0], xKevPM2[1], xKevPM2[2], color = [0.25, 0.25, 0], linewidth = 3)
ax.plot(xKerbin[0], xKerbin[1], xKerbin[2], '--', color = [0, 0, 0.5])
ax.plot(xKevMM[0], xKevMM[1], xKevMM[2], '--', color = [0, 0.5, 0])
plt.show()


fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')

ax.plot(xLr1[0] - xKerbin1[0], xLr1[1] - xKerbin1[1], xLr1[2] - xKerbin1[2], color = [0, 0, 1])
ax.plot(xLr2_[0] - xKerbin2[0], xLr2_[1] - xKerbin2[1], xLr2_[2] - xKerbin2[2], color = [1, 0, 0])
ax.plot(xLr3[0] - xKerbin3[0], xLr3[1] - xKerbin3[1], xLr3[2] - xKerbin3[2], color = [0, 0, 1])
ax.plot(xKev[0] - xKerbin[0], xKev[1] - xKerbin[1], xKev[2] - xKerbin[2], color = [0, 1, 0])
#ax.plot(xKevPM1[0] - xKerbinKev[0], xKevPM[1] - xKerbinKev[1], xKevPM[2] - xKerbinKev[2], color = [0.5, 0.5, 0])
plt.show()

tck.toc("Time : ")