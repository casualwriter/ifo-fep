"""M4 (loop 3 of plan/20261016-k1-math-deep-v1.md): floor law for a broken sensory/active split. Predictions fixed before running:
L3a eps = eta = 0: median conditional transfer ratio ~ delta^p, p (m = 8 -> 16) in [3.6, 4.4], both directions.
L3b floor (m = 32) slope in log eps in [1.8, 2.2] for eps in {0.01, 0.03, 0.1}.
L3c eps = 0, correlated D_B (sensory 2 - active 4): floor slope in log eta in [1.8, 2.2] for eta in {0.01, 0.03, 0.1}.
8 nodes I={0,1}, S={2,3}, A={4,5}, E={6,7}; I hears I,S,A; A hears I,S,A; S hears E,S,A; E hears E,S,A; A[4,0] = 0.5 fixed;
edge prob 0.6, U(0.2,0.8) +-, diag -1, stability shift; seed 9191, 200 nets; x(0) known, push sigma_u = 0.1, h = 2."""
import os, sys, json
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M3_loop2.py"), encoding="utf-8").read().split("\nout = {}")[0])   # load Obs only
I, S, A_, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + A_
hears = {**{k: I + S + A_ for k in I + A_}, **{k: E + S + A_ for k in S + E}}
rng = np.random.default_rng(9191); su = 0.1
def gen():
    M = -np.eye(8)
    for i in range(8):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: M[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    M[4, 0] = 0.5; return M
def stab(M):
    m = np.linalg.eigvals(M).real.max(); return M - (m + 0.1) * np.eye(8) if m > -0.1 else M
def ratio(M, m, D=None):
    o = Obs(M, 2.0, m, D=D); P = su * np.eye(8)
    return o.lc(P[:, I], I, B, E), o.lc(P[:, E], E, B, I)
nets = [gen() for _ in range(200)]; out = {}
# L3a
r = {m: [] for m in [1, 4, 8, 16, 32]}
for M0 in nets:
    M = stab(M0.copy())
    for m in r: r[m].append(ratio(M, m))
med = {m: np.nanmedian(np.array(v), axis=0).tolist() for m, v in r.items()}
p = [float(np.log(med[8][d] / med[16][d]) / np.log(2)) for d in range(2)]
out['L3a'] = dict(median_IE_EI_by_m=med, slope_p_m8_m16=p); print('L3a', out['L3a'], flush=True)
# L3b, L3c
def floors(kind, vals):
    f = {v: [] for v in vals}
    for M0 in nets:
        for v in vals:
            M = M0.copy(); D = None
            if kind == 'eps': M[4, 6] = v
            else: D = np.eye(8); D[2, 4] = D[4, 2] = v
            f[v].append(np.nanmax(ratio(stab(M), 32, D)))
    meds = {v: float(np.nanmedian(x)) for v, x in f.items()}; vs = sorted(vals)
    return dict(median_floor=meds, slope=float(np.polyfit(np.log(vs), np.log([meds[v] for v in vs]), 1)[0]),
                eps0_baseline=float(np.nanmedian([np.nanmax(x) for x in r[32]])))
out['L3b'] = floors('eps', [0.01, 0.03, 0.1]); print('L3b', out['L3b'], flush=True)
out['L3c'] = floors('eta', [0.01, 0.03, 0.1]); print('L3c', out['L3c'], flush=True)
json.dump(out, open(os.path.join(HERE, 'M4.json'), 'w'), indent=1, default=float)
