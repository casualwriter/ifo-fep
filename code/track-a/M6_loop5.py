"""M6 (loop 5 of plan/20261016-k1-math-deep-v1.md): blanket near a bifurcation. Predictions fixed before running:
P1 median slowest-mode share on I at gamma_c >= 0.8.
P2 median total-effect rho_3 (planted, B = S u A) falls by >= 0.2 from gamma = 0 to gamma_c; null (negative 2-cycle) change < 0.05.
P3 median conditional rho_c (h = 3, delta = 0.1, x(0) known, push 0.1) at gamma_c >= 0.95.
Min n: >= 50 nets with real lambda_1 at gamma_c, else inconclusive.
Nets: Loop 3 S/A split, no violations; seed 9393; 100 nets. Internal loop gamma (e0 e1' + e1 e0'); gamma_c by bisection: lambda_1 = -0.02."""
import os, sys, json
import numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M3_loop2.py"), encoding="utf-8").read().split("\nout = {}")[0])
I, S, A_, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + A_
hears = {**{k: I + S + A_ for k in I + A_}, **{k: E + S + A_ for k in S + E}}
rng = np.random.default_rng(9393); su = 0.1
def gen():
    M = -np.eye(8)
    for i in range(8):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: M[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    m = np.linalg.eigvals(M).real.max(); return M - (m + 0.3) * np.eye(8) if m > -0.3 else M
L = np.zeros((8, 8)); L[0, 1] = L[1, 0] = 1.0
K = np.zeros((8, 8)); K[0, 1] = 1.0; K[1, 0] = -1.0   # FIX after run 1: true negative-feedback null (antisymmetric); -L was still a positive cycle
def ldir(R2, src, mid, tgt):
    s = lambda t: R2[np.ix_(t, src)].sum(); d = s(mid) + s(tgt); return np.nan if d < 1e-300 else s(tgt) / d
def rho_tot(M, h): R2 = expm(M * h) ** 2; return 1 - np.nanmax([ldir(R2, I, B, E), ldir(R2, E, B, I)])   # FIX after run 1: 0/0 -> NaN
def slow(M):
    lam, V = np.linalg.eig(M); k = np.argmax(lam.real); v2 = np.abs(V[:, k]) ** 2
    return lam[k], v2
def rho_inf(v2):
    vb, vi, ve = v2[B].sum(), v2[I].sum(), v2[E].sum(); return vb / (vb + max(vi, ve))
def rho_cond(M):
    o = Obs(M, 3.0, 30); P = su * np.eye(8)
    return 1 - np.nanmax([o.lc(P[:, I], I, B, E), o.lc(P[:, E], E, B, I)])
def gamma_c(M0, sign):
    lo, hi = 0.0, 0.05
    while slow(M0 + sign * hi * L)[0].real < -0.02 and hi < 50: hi *= 2
    if hi >= 50: return None
    for _ in range(50):
        mid = (lo + hi) / 2
        if slow(M0 + sign * mid * L)[0].real < -0.02: lo = mid
        else: hi = mid
    return lo
rows = []
for g in range(100):
    M0 = gen(); row = {}
    gcp = gamma_c(M0, 1)
    for Lk, key in [(L, 'pos'), (K, 'neg')]:
        gc = gcp                                        # null uses the same gamma values as the positive loop
        if gc is None: row[key] = None; continue
        r = dict(gc=gc)
        for f in [0, 0.5, 0.9, 1.0]:
            M = M0 + f * gc * Lk; lam, v2 = slow(M)
            r[f] = dict(lam_re=float(lam.real), real=bool(abs(lam.imag) < 1e-9), shareI=float(v2[I].sum() / v2.sum()),
                        rho3=float(rho_tot(M, 3.0)), rho_inf=float(rho_inf(v2)), rho_c=float(rho_cond(M)) if f in (0, 1.0) else None)
        row[key] = r
    rows.append(row); print(g, {k: (v[1.0] if v else None) for k, v in row.items()}, flush=True)
P = [r['pos'] for r in rows if r['pos'] and r['pos'][1.0]['real']]
N = [r['neg'] for r in rows if r['neg']]
med = lambda L_, f, k: float(np.nanmedian([x[f][k] for x in L_]))
out = dict(n_pos_real=len(P), n_neg=len(N),
           pos={f: {k: med(P, f, k) for k in ['lam_re', 'shareI', 'rho3', 'rho_inf']} for f in [0, 0.5, 0.9, 1.0]},
           neg={f: {k: med(N, f, k) for k in ['lam_re', 'shareI', 'rho3']} for f in [0, 0.5, 0.9, 1.0]},
           rho_c_pos={f: med(P, f, 'rho_c') for f in [0, 1.0]})
out['P1'] = out['pos'][1.0]['shareI']; out['P2_drop'] = out['pos'][0]['rho3'] - out['pos'][1.0]['rho3']
out['P2_null_change'] = abs(out['neg'][0]['rho3'] - out['neg'][1.0]['rho3']); out['P3'] = out['rho_c_pos'][1.0]
print(json.dumps(out, indent=1)); json.dump(dict(summary=out, rows=rows), open(os.path.join(HERE, 'M6.json'), 'w'), indent=1, default=float)
