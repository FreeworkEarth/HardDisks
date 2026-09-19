#!/usr/bin/env python3
"""##CHRIS: per-run frequency by a time-domain damped-cosine fit, replacing the
binned FFT peak. The binned estimator resolves nu only to df/nu = 3.2-3.4% (the
route-A traces), which is 5x the 0.70% spread of sqrt(T_i) and therefore cannot
see the seeding-temperature signal at all. A parametric fit uses every sample and
returns a genuine standard error.

Model (post-transient, mean removed):
    x(t) = A exp(-gamma (t-t0)) cos(2 pi nu (t-t0) + phi) + c
Seeded by linear prediction (Prony): a damped sinusoid satisfies
    x[n] = 2 r cos(w dt) x[n-1] - r^2 x[n-2],
so a least-squares solve for (a,b) gives r = sqrt(-b), w = arccos(a/(2r))/dt.
Then curve_fit refines and gives sigma_nu from the covariance.
Analysis only; reads the accepted traces read-only."""
import csv, math, os, sys, warnings
import numpy as np
from scipy.optimize import curve_fit
warnings.filterwarnings("ignore")
DROP = 0.20            # same transient fraction the pipeline uses

def prony_seed(t, x):
    dt = float(np.median(np.diff(t)))
    if len(x) < 8 or dt <= 0: return None
    A = np.vstack([x[1:-1], x[:-2]]).T
    try: ab, *_ = np.linalg.lstsq(A, x[2:], rcond=None)
    except Exception: return None
    a, b = float(ab[0]), float(ab[1])
    if not (-1.0 < b < 0.0): return None
    r = math.sqrt(-b)
    c = a/(2*r)
    if not (-1.0 < c < 1.0): return None
    w = math.acos(c)/dt
    gam = -math.log(max(r, 1e-300))/dt
    return w/(2*math.pi), max(gam, 0.0), dt

def model(t, A, gam, nu, phi, c):
    return A*np.exp(-gam*t)*np.cos(2*math.pi*nu*t + phi) + c

def fit_trace(path, nu_hint=None):
    t=[]; x=[]
    with open(path) as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            try: t.append(float(row["Time"])); x.append(float(row["Displacement(σ)"]))
            except (KeyError, ValueError): pass
    if len(t) < 50: return None
    t=np.asarray(t); x=np.asarray(x)
    i0=int(DROP*len(t)); t=t[i0:]-t[i0]; x=x[i0:]
    x0=x-x.mean()
    s=prony_seed(t,x0)
    if s is None:
        if nu_hint is None: return None
        nu0, gam0 = nu_hint, 1.0/max(t[-1],1e-9)
    else:
        nu0, gam0, _ = s
        if nu_hint and (nu0<=0 or abs(nu0-nu_hint)/nu_hint > 0.5): nu0 = nu_hint
    A0=float(np.std(x0)*math.sqrt(2)) or 1.0
    try:
        p,cov = curve_fit(model, t, x0, p0=[A0, max(gam0,1e-12), nu0, 0.0, 0.0],
                          maxfev=20000)
    except Exception:
        return None
    nu=abs(float(p[2]))
    if not np.all(np.isfinite(cov)): return None
    snu=float(math.sqrt(abs(cov[2,2])))
    resid=x0-model(t,*p)
    return dict(nu=nu, sigma_nu=snu, gamma=abs(float(p[1])),
                rms=float(np.sqrt(np.mean(resid**2))), amp=abs(float(p[0])), n=len(t))

if __name__ == "__main__":
    for f in sys.argv[1:]:
        r=fit_trace(f)
        print(os.path.basename(f), r)
