"""M1b (POST-HOC, after M1): two M1 predictions missed their range.
C3 at K = 4000: no-sham false rejection 0.825 vs formula 0.736 (range 0.08). Check: is it sampling error in b (200k trials) and in the
  200 repeats? Here: b from 2M trials, 1000 repeats at K = 4000, same data (seed 7070).
C5 at s = 0.05: pred / exact = 1.23 (range 10%). Check: is it the O(s) correction (s ||Q|| is not small for Q ~ N(0,1))?
  Here: pred / exact at s = 0.005 and 0.0005 (should -> 1 if the leading term is right)."""
import os, sys, json
import numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
rng = np.random.default_rng(7070); dt, h1 = 0.01, 0.05
A = np.array([[-1.0, 0, 0.8], [0, -1.0, 0.8], [0, 0, -0.1]])
x = np.zeros(3); X = np.empty((10000, 2))
for _ in range(3000): x = x + A @ x * dt + rng.standard_normal(3) * np.sqrt(dt)
for k in range(10000):
    for _ in range(5): x = x + A @ x * dt + rng.standard_normal(3) * np.sqrt(dt)
    X[k] = x[:2]
F1 = np.linalg.lstsq(X[:-1], X[1:], rcond=None)[0]; F2 = np.linalg.lstsq(X[:-2], X[2:], rcond=None)[0]
def trials(K, u):
    idx = rng.integers(0, len(X) - 2, K); s = np.zeros((K, 3)); s[:, :2] = X[idx]; s[:, 2] = rng.standard_normal(K) * np.sqrt(5.0)
    x0 = s[:, :2].copy(); s[:, 0] += u
    for st in range(10):
        s = s + s @ A.T * dt + rng.standard_normal((K, 3)) * np.sqrt(dt)
        if st == 4: r1 = s[:, :2] - x0 @ F1
    return ((4 * r1 - (s[:, :2] - x0 @ F2)) / (2 * h1))[:, 1]
big = np.concatenate([trials(200000, 0.0) for _ in range(10)]); b, sig = big.mean(), big.std()
K, zc, u = 4000, norm.ppf(0.975), 3.0
ns = []
for rep in range(1000):
    p = trials(K, u) / u; ns.append(abs(p.mean() / (p.std(ddof=1) / np.sqrt(K))) > zc)
m = b * np.sqrt(K) / sig; m_se = (sig / np.sqrt(len(big))) * np.sqrt(K) / sig
f = lambda m: norm.sf(zc - m) + norm.cdf(-zc - m)
out['C3'] = dict(b=float(b), sd=float(sig), empirical=float(np.mean(ns)), empirical_se=float(np.std(ns) / np.sqrt(1000)),
                 formula=float(f(m)), formula_range_from_b_se=[float(f(m - 2 * m_se)), float(f(m + 2 * m_se))])
rng = np.random.default_rng(8080); I, B, E = [0, 1], [2, 3], [4, 5]; rat = {0.05: [], 0.005: [], 0.0005: []}
def leak_dir(R2, src, mid, tgt):
    s = lambda t: R2[np.ix_(t, src)].sum(); return s(tgt) / (s(mid) + s(tgt))
for g in range(200):
    M = rng.uniform(-1, 1, (6, 6)) * (rng.random((6, 6)) < 0.6); M = (M + M.T) / 2
    M[np.ix_(I, E)] = 0; M[np.ix_(E, I)] = 0; np.fill_diagonal(M, 0)
    Pi = M + (0.5 - np.linalg.eigvalsh(M).min()) * np.eye(6)
    Q = rng.standard_normal((6, 6)); Q = Q - Q.T; QP = Q @ Pi
    if (Pi[np.ix_(B, I)] ** 2).sum() == 0 or (Pi[np.ix_(B, E)] ** 2).sum() == 0: continue
    c = max((QP[np.ix_(E, I)] ** 2).sum() / (Pi[np.ix_(B, I)] ** 2).sum(), (QP[np.ix_(I, E)] ** 2).sum() / (Pi[np.ix_(B, E)] ** 2).sum())
    for s in rat:
        R2 = (-(np.eye(6) + s * Q) @ Pi * (1 - np.eye(6))) ** 2
        rat[s].append(s ** 2 * c / max(leak_dir(R2, I, B, E), leak_dir(R2, E, B, I)))
out['C5'] = {str(s): float(np.median(v)) for s, v in rat.items()}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, 'M1b.json'), 'w'), indent=1)
