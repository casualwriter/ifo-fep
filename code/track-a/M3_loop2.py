"""M3 (loop 2 of plan/20261016-k1-math-deep-v1.md): conditional, unit-invariant crossing share. Predictions fixed before running.
Observer: knows x(0), pushes I with u ~ N(0, I), reads B at t_k = k h / m (k = 1..m) and E at h.
  lc_IE = I(u; x_E(h) | x(0), B samples) / I(u; x_E(h) | x(0)); rho_c = 1 - max(lc_IE, lc_EI). Exact Gaussian algebra:
  x(t) - e^{At} x(0) = e^{At} P u + N(t), Cov N(t_k), N(t_l) = e^{A(t_k - t_l)} W(t_l), W(t) = Sig - e^{At} Sig e^{A't}, A Sig + Sig A' + D = 0.
  I(u; S) = 1/2 [logdet(M_S M_S' + CN_S) - logdet(CN_S)].
L2a: X2 wiring cuts (seed 21, 200 nets), h = 2, m = 1, 2, 4, 8, 16. Prediction: median lc falls monotonically in m and
     median(m = 16) < 0.1 * median(m = 1).
L2b: T4 nets (seed 5050, first 50), direct-coupling rho (O1). S = diag(10^U(0,1)). Prediction: best candidate changes in > 10% of nets.
     Invariance: rho_c with (S A S^-1, S D S, push S_I) equals rho_c of the original to < 1e-9 (10 candidates x 5 nets).
L2c (main): first 20 T4 nets, 2000 random candidates + O1's best. Fast (h = 0.3, m = 4) vs slow (h = 3, m = 4).
     Prediction: best rho_c candidate differs in >= 30% of nets. Death clause: same in >= 90% -> h dependence was a total-effect artefact."""
import os, sys, json
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'causal-partition-v0_6'))
import expC
from T1_graded_blankets import n as N12, net, PARTS
from T4_observer_dependence import rho_from
ld = lambda M: np.linalg.slogdet(M)[1]
class Obs:
    def __init__(self, A, h, m, D=None):
        k = len(A); D = np.eye(k) if D is None else D; self.k, self.m = k, m
        Sig = solve_continuous_lyapunov(A, -D); ts = h * np.arange(1, m + 1) / m
        self.Phi = [expm(A * t) for t in ts]
        W = [Sig - P @ Sig @ P.T for P in self.Phi]
        self.CN = np.zeros((k * m, k * m))
        for a in range(m):
            for b in range(a + 1):
                blk = (expm(A * (ts[a] - ts[b])) @ W[b])
                self.CN[a * k:(a + 1) * k, b * k:(b + 1) * k] = blk; self.CN[b * k:(b + 1) * k, a * k:(a + 1) * k] = blk.T
    def lc(self, P, src_cols, mid, tgt):
        """P: k x p push matrix (columns = pushed directions). mid, tgt: node lists. Returns conditional share."""
        k, m = self.k, self.m
        M = np.vstack([Ph @ P for Ph in self.Phi])
        z = [a * k + b for a in range(m) for b in mid]; y = [(m - 1) * k + e for e in tgt]
        I_ = lambda S: 0.5 * (ld(M[S] @ M[S].T + self.CN[np.ix_(S, S)]) - ld(self.CN[np.ix_(S, S)]))
        tot = I_(y)
        return np.nan if tot < 1e-12 else (I_(y + z) - I_(z)) / tot
    def rho_c(self, I, B, E, S=None):
        S = np.eye(self.k) if S is None else S
        v = [self.lc(S[:, I], I, B, E), self.lc(S[:, E], E, B, I)]
        return np.nan if np.isnan(v).all() else 1 - np.nanmax(v)
out = {}

# L2a
expC.rng = np.random.default_rng(21); L = {m: [] for m in [1, 2, 4, 8, 16]}
for g in range(200):
    A, I, B, E = expC.random_planted(0.0); A6 = A[:6, :6]
    for m in L:
        o = Obs(A6, 2.0, m); v = [o.lc(np.eye(6)[:, I], I, B, E), o.lc(np.eye(6)[:, E], E, B, I)]
        L[m].append(np.nanmax(v) if not np.isnan(v).all() else np.nan)
med = {m: float(np.nanmedian(v)) for m, v in L.items()}
out['L2a'] = dict(median_lc_by_m=med, q90_by_m={m: float(np.nanquantile(v, 0.9)) for m, v in L.items()},
                  loglog_slope_m8_m16=float(np.log(med[16] / med[8]) / np.log(0.5)))
print('L2a', out['L2a'], flush=True)

# L2b
rng = np.random.default_rng(5050); nets = []
for g in range(50):
    A = net(rng, 0.0); rng.choice(N12, 2, replace=False); nets.append(A[:N12, :N12].copy())
srng = np.random.default_rng(6060); flips = []; sp = []
for Ao in nets:
    W = Ao * (1 - np.eye(N12)); s = 10 ** srng.uniform(0, 1, N12)
    r0 = rho_from(W ** 2, range(N12)); r1 = rho_from((s[:, None] * W / s[None, :]) ** 2, range(N12))
    ok = ~np.isnan(r0) & ~np.isnan(r1)
    flips.append(int(np.argmax(np.where(ok, r0, -1)) != np.argmax(np.where(ok, r1, -1)))); sp.append(float(spearmanr(r0[ok], r1[ok]).statistic))
inv = []
for Ao in nets[:5]:
    s = 10 ** srng.uniform(0, 1, N12); S = np.diag(s); Si = np.diag(1 / s)
    o0 = Obs(Ao, 1.0, 2); o1 = Obs(S @ Ao @ Si, 1.0, 2, D=S @ S)
    for p in srng.choice(len(PARTS), 10, replace=False):
        I, E, B = [list(x) for x in PARTS[p]]
        inv.append(abs(o0.rho_c(I, B, E) - o1.rho_c(I, B, E, S=S)))
out['L2b'] = dict(old_rho_best_flips=float(np.mean(flips)), old_rho_spearman_median=float(np.median(sp)), rho_c_max_abs_change=float(np.nanmax(inv)))
print('L2b', out['L2b'], flush=True)

# L2c
crng = np.random.default_rng(7171); rows = []
for g, Ao in enumerate(nets[:20]):
    W = Ao * (1 - np.eye(N12)); r0 = rho_from(W ** 2, range(N12))
    cand = list(crng.choice(len(PARTS), 2000, replace=False)) + [int(np.nanargmax(r0))]
    of, os_ = Obs(Ao, 0.3, 4), Obs(Ao, 3.0, 4)
    rf, rs = [], []
    for p in cand:
        I, E, B = [list(x) for x in PARTS[p]]; rf.append(of.rho_c(I, B, E)); rs.append(os_.rho_c(I, B, E))
    rf, rs = np.array(rf), np.array(rs); ok = ~np.isnan(rf) & ~np.isnan(rs)
    bf, bs = np.argmax(np.where(ok, rf, -9)), np.argmax(np.where(ok, rs, -9))
    k = 20; tf, ts = set(np.argsort(-np.where(ok, rf, -9))[:k]), set(np.argsort(-np.where(ok, rs, -9))[:k])
    rows.append(dict(differs=bool(bf != bs), spearman=float(spearmanr(rf[ok], rs[ok]).statistic), top1pct_overlap=len(tf & ts) / k,
                     best_fast=float(rf[bf]), slow_of_fastbest=float(rs[bf]), best_slow=float(rs[bs]), spearman_fast_vs_r0=float(spearmanr(rf[ok], r0[cand][ok]).statistic)))
    print(g, rows[-1], flush=True)
out['L2c'] = dict(differs=float(np.mean([r['differs'] for r in rows])), **{k: float(np.median([r[k] for r in rows])) for k in rows[0] if k != 'differs'}, rows=rows)
print('L2c', {k: v for k, v in out['L2c'].items() if k != 'rows'})
json.dump(out, open(os.path.join(HERE, 'M3.json'), 'w'), indent=1, default=float)
