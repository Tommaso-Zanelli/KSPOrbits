#!/usr/bin/env python3
import numpy as np

def ofTest(t1, t2, dv0, verbose = False, maxIter = 2147483648):
    return 1.0 - 0.5*t1 + 1.5*t2 + 5*t1*t1 - 0.2*t1*t2 + 4*t2*t2, np.zeros(3), 0


def optimizeCost(of, t1, t2, dvn, verbose = False, maxIter = 2147483648):
    # Setup parameters
    minRelCGain = np.exp2(-26) # Fraction of initial cost below which iteration is stopped 
    DMax = 100                 # Maximum distance from target allowable for convergence
    largeRel = 1048576.0       # Ratio of derivatives in t1 and t2 above which one of the two is ignored
    tmin = 1                   # Minimum reference variable for computing the steps
    dt2N = 21                  # Power of two ratio between variable and step
    dtMin2N = 48               # Power of two ratio between variable and minimum step
    aMin = np.exp2(-26)        # Minimum increase/decrease along gradient direction
    maxIter = 1024             # Hardcoded maximum number of iterations

    outerIterate = True
    Ci = None
    nIter = 0
    while outerIterate:
        # Reference values at current position
        C0, dv0, D0 = of(t1, t2, dv0 = dvn, verbose = verbose, maxIter = maxIter)
        if Ci == None:
            Ci = C0
            Cf = Ci

        # 2-base exponent of the input (capped)
        t1ev = np.floor(np.log2(max(abs(t1), tmin)))
        t2ev = np.floor(np.log2(max(abs(t2), tmin)))

        # Time steps for computing finite difference derivatives
        dt1 = np.exp2(t1ev - dt2N)
        dt2 = np.exp2(t2ev - dt2N)
        
        # Time steps inverses
        dt1m1 = 1.0/dt1
        dt2m1 = 1.0/dt2

        # Minimum time step sizes
        dt1Min = np.exp2(t1ev - dtMin2N)
        dt2Min = np.exp2(t2ev - dtMin2N)

        # Compute values at increased and decreased times
        if verbose: print(f"Computing derivatives for iteration {nIter:d}")
        iterate       = True
        p1converged   = False
        m1converged   = False
        p2converged   = False
        m2converged   = False
        p1p2converged = False
        m1m2converged = False
        m1p2converged = False
        p1m2converged = False
        while iterate:
            if not p1converged:
                Cp1, dvp1, Dp1 = of(t1 + dt1, t2, dv0 = dvn, verbose = verbose, maxIter = maxIter)
                p1converged = (Dp1 < DMax and not np.isnan(Dp1))
            if (not m1converged) and p1converged:
                Cm1, dvm1, Dm1 = of(t1 - dt1, t2, dv0 = (2.0*dvn - dvp1), verbose = verbose, maxIter = maxIter)
                m1converged = (Dm1 < DMax and not np.isnan(Dm1))
            if not p2converged:
                Cp2, dvp2, Dp2 = of(t1, t2 + dt2, dv0 = dvn, verbose = verbose, maxIter = maxIter)
                p2converged = (Dp2 < DMax and not np.isnan(Dp2))
            if (not m2converged) and p2converged:
                Cm2, dvm2, Dm2 = of(t1, t2 - dt2, dv0 = (2.0*dvn - dvp2), verbose = verbose, maxIter = maxIter)
                m2converged = (Dm2 < DMax and not np.isnan(Dm2))
            if (not p1p2converged) and p1converged and p2converged:
                Cp1p2, dvp1p2, Dp1p2 = of(t1 + dt1, t2 + dt2, dv0 = (dvp1 + dvp2 - dvn), verbose = verbose, maxIter = maxIter)
                p1p2converged = (Dp1p2 < DMax and not np.isnan(Dp1p2))
            if (not m1m2converged) and p1converged and p2converged and p1p2converged:
                Cm1m2, dvm1m2, Dm1m2 = of(t1 - dt1, t2 - dt2, dv0 = (dvm1 + dvm2 - dvn), verbose = verbose, maxIter = maxIter)
                m1m2converged = (Dm1m2 < DMax and not np.isnan(Dm1m2))
            if (not p1m2converged) and p1converged and p2converged and p1p2converged and m1m2converged:
                Cp1m2, dvp1m2, Dp1m2 = of(t1 + dt1, t2 - dt2, dv0 = (dvp1 + dvm2 - dvn), verbose = verbose, maxIter = maxIter)
                p1m2converged = (Dp1m2 < DMax and not np.isnan(Dp1m2))
            if (not m1p2converged) and p1converged and p2converged and p1p2converged and m1m2converged and p1m2converged:
                Cm1p2, dvm1p2, Dm1p2 = of(t1 - dt1, t2 + dt2, dv0 = (dvm1 + dvp2 - dvn), verbose = verbose, maxIter = maxIter)
                m1p2converged = (Dm1p2 < DMax and not np.isnan(Dm1p2))
            allMixedConverged = p1p2converged and m1m2converged and p1m2converged and m1p2converged
            t1OnlyConverged = p1converged and m1converged
            t2OnlyConverged = p2converged and m2converged
            allConverged = t1OnlyConverged and allMixedConverged and t2OnlyConverged
            if (not t1OnlyConverged) or (not allMixedConverged):
                dt1 *= 0.5
                p1converged   = False
                m1converged   = False
                p1p2converged = False
                m1m2converged = False
                p1m2converged = False
                m1p2converged = False
            if (not t2OnlyConverged) or (not allMixedConverged):
                dt2 *= 0.5
                p2converged   = False
                m2converged   = False
                p1p2converged = False
                m1m2converged = False
                p1m2converged = False
                m1p2converged = False
            if allConverged or dt1 < dt1Min or dt2 < dt2Min:
                iterate = False

        # Catch failed convergence while computing increased/decreased values
        if (not allConverged) and (dt1 < dt1Min or dt2 < dt2Min):
            outerIterate = False
            if verbose: print(f"\nStopping at iteration {nIter:d} due to failed convergence while computing derivatives\n")
        else:
            # Compute cost function first and second order derivatives
            dCdt1     = 0.5*(Cp1 - Cm1)*dt1m1
            dCdt2     = 0.5*(Cp2 - Cm2)*dt2m1
            d2Cdt12   = (Cp1 - 2.0*C0 + Cm1)*dt1m1*dt1m1
            d2Cdt22   = (Cp2 - 2.0*C0 + Cm2)*dt2m1*dt2m1
            d2Cdt1dt2 = 0.25*(Cp1p2 - Cp1m2 - Cm1p2 + Cm1m2)*dt1m1*dt2m1
    
            # Compute velocity first and second order derivatives
            ddvdt1     = 0.5*(dvp1 - dvm1)*dt1m1
            ddvdt2     = 0.5*(dvp2 - dvm2)*dt2m1
            d2dvdt12   = (dvp1 - 2.0*dv0 + dvm1)*dt1m1*dt1m1
            d2dvdt22   = (dvp2 - 2.0*dv0 + dvm2)*dt2m1*dt2m1
            d2dvdt1dt2 = 0.25*(dvp1p2 - dvp1m2 - dvm1p2 + dvm1m2)*dt1m1*dt2m1
    
            # Compute time increases for reducing the cost function
            if abs(dCdt1) > abs(dCdt2)*largeRel and abs(d2Cdt12) > abs(d2Cdt22)*largeRel:
                if verbose: print(f"\tVariation along t2 comparatively small (dCdt2/dCdt1 = {dCdt2/dCdt1:g}, d2Cdt22/d2Cdt12 = {d2Cdt22/d2Cdt12:g}), optimizing along t1 only")
                ddt1 = -dCdt1/d2Cdt12
                ddt2 = 0.0
            elif abs(dCdt1)*largeRel < abs(dCdt2) and abs(d2Cdt12)*largeRel < abs(d2Cdt22):
                if verbose: print(f"\tVariation along t1 comparatively small (dCdt1/dCdt2 = {dCdt1/dCdt2:g}, d2Cdt12/d2Cdt22 = {d2Cdt12/d2Cdt22:g}), optimizing along t2 only")
                ddt1 = 0.0
                ddt2 = -dCdt2/d2Cdt22
            else:
                DM = d2Cdt12*d2Cdt22 - d2Cdt1dt2**2
                if verbose: print(f"\tCost Jacobian matrix determinant: {DM:g}")
                if DM == 0:
                    ddt1 = 0.0
                    ddt2 = 0.0
                    outerIterate = False
                    if verbose: print(f"\nStopping at iteration {nIter:d} due to singular Jacobian matrix\n")
                else:
                    DMm1 = 1.0/DM
                    ddt1 = (d2Cdt1dt2*dCdt2 - d2Cdt22*dCdt1)*DMm1
                    ddt2 = (d2Cdt1dt2*dCdt1 - d2Cdt12*dCdt2)*DMm1
    
            # Stop iterating if time increases 
            if verbose: print(f"Finding optimum for iteration {nIter:d}")
            if abs(ddt1) < dt1Min and abs(ddt2) < dt2Min:
                if verbose: print(f"\nStopping at iteration {nIter:d} due to time increases below treshold\n")
                outerIterate = False
            else:
                a = 1.0
                iterate = True
                while iterate:
                    if verbose: print(f"\tTesting amplitude value {a:g}")
                    addt1 = a*ddt1
                    addt2 = a*ddt2
                    ddv = dvn + ddvdt1*addt1 + ddvdt2*addt2 + \
                          0.5*d2dvdt12*addt1*addt1 + d2dvdt1dt2*addt1*addt2 + \
                          0.5*d2dvdt22*addt2*addt2 
                    Cgn, dvgn, Dgn = of(t1 + addt1, t2 + addt2, dv0 = ddv, verbose = verbose, maxIter = maxIter)
                    dC = C0 - Cgn
                    betterResult = (dC > 0 and Dgn < DMax and not np.isnan(Dgn))
                    if betterResult or a < aMin:
                        iterate = False
                    else:
                        a *= 0.5
    
                if betterResult:
                    Cf = Cgn
                    t1 += addt1
                    t2 += addt2
                    dvn = dvgn
                    if verbose: print(f"After iteration {nIter:d} cost went from {C0:g} to {Cf:g}, {100*(C0-Cf)/C0:8.4g}% reduction\n")
                else:
                    outerIterate = False
                    if verbose: print(f"Stopping at iteration {nIter:d} due to failure to decrease cost\n")
                if not abs(dC) > minRelCGain*C0:
                    outerIterate = False
                    if verbose: print(f"Stopping at iteration {nIter:d} due to cost decrease below treshold\n")

        nIter += 1
        if nIter > maxIter:
            outerIterate = False
            if verbose: print(f"Stopping at iteration {nIter:d} as the maximum number of iterations allowed was reached\n")

    return t1, t2, dvn, Ci, Cf

def __main__():
    print(optimizeCost(ofTest, 1, 1, np.zeros(3), verbose = True))
