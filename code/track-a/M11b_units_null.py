"""M11b: discriminating units check (replaces the M11 inference). Prospectively specified.
First 50 T4 networks (seed 5050 stream); for each, 3 random diagonal unit systems S (log-uniform, one decade, seed 6262).
Within each unit system, total-effect share: best partition at h = 0.5 vs h = 3 (main), and null flips h vs 1.2h (0.5 vs 0.6, 3 vs 3.6).
Prediction: main flip rate >= 3x the mean null flip rate, in random units and in noise units.
Death: main <= 1.5x null in random units -> the time-scale effect of the total-effect score is not separable from unit effects."""
import os, sys, json
import numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from T1_graded_blankets import n, net
from T4_observer_dependence import rho_from
rng = np.random.default_rng(5050); srng = np.random.default_rng(6262)
HS = [0.5, 0.6, 3.0, 3.6]; res = {'noise': [], 'random': []}
for g in range(50):
    A = net(rng, 0.0); rng.choice(n, 2, replace=False); Ao = A[:n, :n]
    R = {h: expm(h * Ao) for h in HS}
    systems = [('noise', np.ones(n))] + [('random', 10 ** srng.uniform(0, 1, n)) for _ in range(3)]
    for key, s in systems:
        best = {h: int(np.nanargmax(rho_from((s[:, None] * R[h] / s[None, :]) ** 2, range(n)))) for h in HS}
        res[key].append(dict(main=best[0.5] != best[3.0], null=[best[0.5] != best[0.6], best[3.0] != best[3.6]]))
out = {}
for k, v in res.items():
    m = float(np.mean([x['main'] for x in v])); nl = float(np.mean([y for x in v for y in x['null']]))
    out[k] = dict(n=len(v), main_flip=m, null_flip=nl, ratio=(m / nl if nl else None))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M11b.json'), 'w'), indent=1)
