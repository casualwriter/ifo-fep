"""M5 (loop 4 of plan/20261016-k1-math-deep-v1.md): fair observer-dependence test + minimax objectivity. Predictions fixed before running:
L4a at fixed delta, best eligible candidate differs between h = 0.5 and h = 3 in >= 30% of nets AND >= 3x the null flip rate
    (h vs 1.1 h, mean over h in {0.5, 3}). Death: <= 1.5x null.
L4b flip rate (0.5 vs 3) at delta = 0.1 < at delta = 0.25.
L4c minimax-best candidate (max over candidates of min rho_c over the 6 observers h in {0.5,1,3} x delta in {0.1,0.25}) is the planted
    partition in >= 70% of nets.
Eligible = rho_c >= 0.7 under both observers compared; a net counts if >= 2 eligible.
Nets: Loop 3 S/A generator (I={0,1}, S={2,3}, A={4,5}, E={6,7}) + each forbidden edge with prob 0.1, weight U(0.05,0.2) +-;
40 nets, seed 9292; 210 candidates (|I| = |E| = 2, |B| = 4, I/E unordered). Observer: x(0) known, push 0.1, m = round(h / delta)."""
import os, sys, json, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M3_loop2.py"), encoding="utf-8").read().split("\nout = {}")[0])   # load Obs only
I0, S0, A0, E0 = [0, 1], [2, 3], [4, 5], [6, 7]
hears = {**{k: I0 + S0 + A0 for k in I0 + A0}, **{k: E0 + S0 + A0 for k in S0 + E0}}
rng = np.random.default_rng(9292); su = 0.1
def gen():
    M = -np.eye(8)
    for i in range(8):
        for j in range(8):
            if i == j: continue
            if j in hears[i]:
                if rng.random() < 0.6: M[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
            elif rng.random() < 0.1: M[i, j] = rng.uniform(0.05, 0.2) * rng.choice([-1, 1])
    m = np.linalg.eigvals(M).real.max()
    return M - (m + 0.1) * np.eye(8) if m > -0.1 else M
CANDS = []
for I in itertools.combinations(range(8), 2):
    for E in itertools.combinations([k for k in range(8) if k not in I], 2):
        if I < E: CANDS.append((list(I), [k for k in range(8) if k not in I + E], list(E)))
PLANT = next(c for c, (I, B, E) in enumerate(CANDS) if I == I0 and E == E0)
OBS = [(h, d) for d in [0.1, 0.25] for h in [0.5, 0.55, 1, 3, 3.3]]
def rho_c(o, I, B, E):
    P = su * np.eye(8); v = [o.lc(P[:, I], I, B, E), o.lc(P[:, E], E, B, I)]
    return np.nan if np.isnan(v).all() else 1 - np.nanmax(v)
rows = []
for g in range(40):
    A = gen(); R = {}
    for h, d in OBS:
        o = Obs(A, h, max(1, round(h / d))); R[(h, d)] = np.array([rho_c(o, *c) for c in CANDS])
    def best(o1, o2):
        ok = (R[o1] >= 0.7) & (R[o2] >= 0.7)
        if ok.sum() < 2: return None
        return bool(np.argmax(np.where(ok, R[o1], -9)) != np.argmax(np.where(ok, R[o2], -9)))
    row = {}
    for d in [0.1, 0.25]:
        row[f'flip_{d}'] = best((0.5, d), (3, d)); row[f'null_{d}'] = [best((0.5, d), (0.55, d)), best((3, d), (3.3, d))]
    six = np.min(np.vstack([R[(h, d)] for h, d in OBS if h in (0.5, 1, 3)]), axis=0)
    row.update(minimax_is_plant=bool(np.nanargmax(six) == PLANT), n_pass_all=int(np.nansum(six >= 0.7)), plant_min=float(six[PLANT]),
               minimax_value=float(np.nanmax(six)), plant_rho={f'{h}_{d}': float(R[(h, d)][PLANT]) for h, d in OBS})
    rows.append(row); print(g, {k: v for k, v in row.items() if k != 'plant_rho'}, flush=True)
def rate(key, idx=None):
    v = [r[key] if idx is None else r[key][idx] for r in rows]; v = [x for x in v if x is not None]
    return (float(np.mean(v)) if v else None, len(v))
out = {}
for d in [0.1, 0.25]:
    f, nf = rate(f'flip_{d}'); n1, _ = rate(f'null_{d}', 0); n2, _ = rate(f'null_{d}', 1)
    nul = np.nanmean([x for x in [n1, n2] if x is not None])
    out[f'delta_{d}'] = dict(flip=f, n=nf, null_h05=n1, null_h3=n2, null_mean=float(nul), ratio=(f / nul) if nul else None)
out['L4c'] = dict(minimax_is_plant=float(np.mean([r['minimax_is_plant'] for r in rows])), n_pass_all_median=float(np.median([r['n_pass_all'] for r in rows])),
                  plant_min_median=float(np.median([r['plant_min'] for r in rows])), minimax_value_median=float(np.median([r['minimax_value'] for r in rows])))
print(json.dumps(out, indent=1))
json.dump(dict(summary=out, rows=rows), open(os.path.join(HERE, 'M5.json'), 'w'), indent=1, default=float)
