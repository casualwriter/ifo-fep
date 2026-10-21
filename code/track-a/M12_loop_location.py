"""M12: matched-slowing control for critical slowing (paper Sec. 5.4). Prospectively specified.
Same networks as M6 (seed 9393 generator, 100 nets). Positive 2-cycle placed (a) inside I (nodes 0,1: as M6) or (b) inside the
blanket's sensory pair (nodes 2,3), each raised by bisection until the slowest eigenvalue is -0.02 (same slowness). M stays 0 in both.
Prediction (Prop 4): total-effect rho(3) falls for (a) and does NOT fall for (b): median change (b) >= -0.02 (no drop), while (a)
drops >= 0.2 (replicates M6). Slow-mode share on B at gamma_c for (b) >= 0.6."""
import os, json
import numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
I, S, A_, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + A_
hears = {**{k: I + S + A_ for k in I + A_}, **{k: E + S + A_ for k in S + E}}
rng = np.random.default_rng(9393)
def gen():
    M = -np.eye(8)
    for i in range(8):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: M[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    m = np.linalg.eigvals(M).real.max(); return M - (m + 0.3) * np.eye(8) if m > -0.3 else M
def ldir(R2, src, mid, tgt):
    s_ = lambda t: R2[np.ix_(t, src)].sum(); d = s_(mid) + s_(tgt); return np.nan if d < 1e-300 else s_(tgt) / d
def rho3(M): R2 = expm(3 * M) ** 2; return 1 - np.nanmax([ldir(R2, I, B, E), ldir(R2, E, B, I)])
def slow(M):
    lam, V = np.linalg.eig(M); k = np.argmax(lam.real); return lam[k], np.abs(V[:, k]) ** 2
def gc(M0, L):
    lo, hi = 0.0, 0.05
    while slow(M0 + hi * L)[0].real < -0.02 and hi < 50: hi *= 2
    if hi >= 50: return None
    for _ in range(50):
        mid = (lo + hi) / 2
        if slow(M0 + mid * L)[0].real < -0.02: lo = mid
        else: hi = mid
    return lo
LI = np.zeros((8, 8)); LI[0, 1] = LI[1, 0] = 1.0
LB = np.zeros((8, 8)); LB[2, 3] = LB[3, 2] = 1.0
res = {'I': [], 'B': []}
for g in range(100):
    M0 = gen()
    for key, L in [('I', LI), ('B', LB)]:
        c = gc(M0, L)
        if c is None: continue
        M1 = M0 + c * L; lam, v2 = slow(M1)
        if abs(lam.imag) > 1e-9: continue
        res[key].append(dict(d_rho3=float(rho3(M1) - rho3(M0)), shareB=float(v2[B].sum() / v2.sum()), shareI=float(v2[I].sum() / v2.sum())))
out = {k: dict(n=len(v), median_d_rho3=float(np.nanmedian([x['d_rho3'] for x in v])), median_shareB=float(np.median([x['shareB'] for x in v])),
               median_shareI=float(np.median([x['shareI'] for x in v]))) for k, v in res.items()}
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M12.json'), 'w'), indent=1)
