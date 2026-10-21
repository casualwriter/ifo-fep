"""T4b (POST-HOC, after T4): is the 99.5% "best blanket differs" result driven by near-ties among ~17k candidates?
Same generator and observers as T4 (first 50 networks of seed 5050). Per network: Spearman correlation of rho across candidates
between O1 (direct) and O3 (h = 3); rho under O3 of O1's best, as a quantile of all candidates under O3; overlap of the top 1%."""
import os, sys, json
import numpy as np
from scipy.linalg import expm
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from T1_graded_blankets import n, net
from T4_observer_dependence import rho_from
rng = np.random.default_rng(5050); rows = []
for g in range(50):
    A = net(rng, 0.0); W = A[:n, :n] * (1 - np.eye(n)); rng.choice(n, 2, replace=False)   # keep the T4 random stream
    r1 = rho_from(W ** 2, range(n)); r3 = rho_from(expm(3 * A)[:n, :n] ** 2, range(n))
    ok = ~(np.isnan(r1) | np.isnan(r3)); r1, r3 = r1[ok], r3[ok]
    b1 = np.argmax(r1); q = float((r3 < r3[b1]).mean())
    k = max(1, len(r1) // 100); t1 = set(np.argsort(-r1)[:k]); t3 = set(np.argsort(-r3)[:k])
    rows.append(dict(spearman=float(spearmanr(r1, r3).statistic), o1best_quantile_under_o3=q,
                     top1pct_overlap=len(t1 & t3) / k, rho1_best=float(r1[b1]), rho3_of_o1best=float(r3[b1]), rho3_best=float(r3.max())))
out = {k: float(np.median([r[k] for r in rows])) for k in rows[0]}
print(json.dumps(out, indent=1))
json.dump(dict(median=out, rows=rows), open(os.path.join(HERE, 'T4b.json'), 'w'), indent=1)
