"""M7: fair test of 'the best partition depends on h' (plan/20261016-k1-next-steps-v1.md section 1). Predictions fixed before running:
F1 best eligible candidate differs between h = 0.5 and h = 3 (delta = 0.1) in >= 30% of eligible nets AND >= 3x the null flip rate
   (mean of 0.5 vs 0.6 and 3 vs 3.6). Death: <= 1.5x null. Min n: >= 30 eligible nets, else inconclusive.
F2 flip rate in 'direct' nets > in 'split' nets.
Eligible: top 10% of rho_c under both observers compared; a net counts if >= 2 eligible.
Nets: S/A generator (I={0,1}, S={2,3}, A={4,5}, E={6,7}), seed 9494, 150 nets, types alternate: 'direct' (I-E edges) or 'split'
(S hears I, A hears E); violation prob 0.25 per pair, weight U(0.05, 0.2) +-. Observer: x(0) known, push 0.1, m = h / 0.1."""
import os, sys, json, itertools
import numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M3_loop2.py"), encoding="utf-8").read().split("\nout = {}")[0])
I0, S0, A0, E0 = [0, 1], [2, 3], [4, 5], [6, 7]
hears = {**{k: I0 + S0 + A0 for k in I0 + A0}, **{k: E0 + S0 + A0 for k in S0 + E0}}
viol = {'direct': [(i, e) for i in I0 for e in E0] + [(e, i) for i in I0 for e in E0],     # (target, source)
        'split': [(s, i) for s in S0 for i in I0] + [(a, e) for a in A0 for e in E0]}
rng = np.random.default_rng(9494); su = 0.1
def gen(kind):
    M = -np.eye(8)
    for i in range(8):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: M[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    for (t, s) in viol[kind]:
        if rng.random() < 0.25: M[t, s] = rng.uniform(0.05, 0.2) * rng.choice([-1, 1])
    m = np.linalg.eigvals(M).real.max(); return M - (m + 0.1) * np.eye(8) if m > -0.1 else M
CANDS = []
for I in itertools.combinations(range(8), 2):
    for E in itertools.combinations([k for k in range(8) if k not in I], 2):
        if I < E: CANDS.append((list(I), [k for k in range(8) if k not in I + E], list(E)))
PLANT = next(c for c, (I, B, E) in enumerate(CANDS) if I == I0 and E == E0)
HS = [0.5, 0.6, 3.0, 3.6]
def rho_c(o, I, B, E):
    P = su * np.eye(8); v = [o.lc(P[:, I], I, B, E), o.lc(P[:, E], E, B, I)]
    return np.nan if np.isnan(v).all() else 1 - np.nanmax(v)
def rho_tot(A, h, I, B, E):
    R2 = expm(A * h) ** 2
    def ld(src, tgt):
        s = lambda t: R2[np.ix_(t, src)].sum(); d = s(B) + s(tgt); return np.nan if d < 1e-300 else s(tgt) / d
    v = [ld(I, E), ld(E, I)]; return np.nan if np.isnan(v).all() else 1 - np.nanmax(v)
def compare(R1, R2):
    ok = ~np.isnan(R1) & ~np.isnan(R2)
    el = ok & (R1 >= np.nanquantile(R1, 0.9)) & (R2 >= np.nanquantile(R2, 0.9))
    if el.sum() < 2: return None
    a, b = np.where(el, R1, -9), np.where(el, R2, -9); ia, ib = np.argmax(a), np.argmax(b)
    sa = np.sort(a[el])[::-1]
    return dict(flip=bool(ia != ib), margin=float(sa[0] - sa[1]))
rows = []
for g in range(150):
    kind = 'direct' if g % 2 == 0 else 'split'; A = gen(kind)
    OB = {h: Obs(A, h, int(round(h / 0.1))) for h in HS}
    Rc = {h: np.array([rho_c(OB[h], *c) for c in CANDS]) for h in HS}
    Rt = {h: np.array([rho_tot(A, h, *c) for c in CANDS]) for h in [0.5, 0.6, 3.0, 3.6]}
    row = dict(kind=kind, main=compare(Rc[0.5], Rc[3.0]), null=[compare(Rc[0.5], Rc[0.6]), compare(Rc[3.0], Rc[3.6])],
               tot_main=compare(Rt[0.5], Rt[3.0]), tot_null=[compare(Rt[0.5], Rt[0.6]), compare(Rt[3.0], Rt[3.6])],
               plant_best={h: bool(np.nanargmax(Rc[h]) == PLANT) for h in HS}, plant_rho={h: float(Rc[h][PLANT]) for h in HS})
    rows.append(row); print(g, kind, row['main'], [x and x['flip'] for x in row['null']], flush=True)
def rate(rs, key, idx=None):
    v = [(r[key] if idx is None else r[key][idx]) for r in rs]; v = [x['flip'] for x in v if x is not None]
    return (float(np.mean(v)) if v else None, len(v))
out = {}
for name, rs in [('all', rows), ('direct', [r for r in rows if r['kind'] == 'direct']), ('split', [r for r in rows if r['kind'] == 'split'])]:
    f, nf = rate(rs, 'main'); n0, _ = rate(rs, 'null', 0); n1, _ = rate(rs, 'null', 1)
    tf, _ = rate(rs, 'tot_main'); t0, _ = rate(rs, 'tot_null', 0); t1, _ = rate(rs, 'tot_null', 1)
    nul = np.mean([x for x in [n0, n1] if x is not None]); tn = np.mean([x for x in [t0, t1] if x is not None])
    out[name] = dict(flip=f, n=nf, null=float(nul), ratio=(f / nul if nul else None), tot_flip=tf, tot_null=float(tn),
                     plant_best={h: float(np.mean([r['plant_best'][h] for r in rs])) for h in HS},
                     margin_median=float(np.median([r['main']['margin'] for r in rs if r['main']])) if nf else None)
print(json.dumps(out, indent=1))
json.dump(dict(summary=out, rows=rows), open(os.path.join(HERE, 'M7.json'), 'w'), indent=1, default=float)
