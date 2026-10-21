"""M11: is the 99.5% 'fast vs slow best partition differs' (T4, total-effect share, noise units) a unit artefact?
Prospectively specified. First 50 T4 networks (seed 5050 stream). For each, one random diagonal S (log-uniform, one decade, seed 6161);
compute the fast share (direct couplings) and the slow share (h = 3) both in the rescaled units (J -> S J S^-1, R -> S R S^-1).
Prediction: fast-best != slow-best in >= 90% of networks in random units (h effect survives units).
Also report: fast-best in noise units vs random units (unit effect alone)."""
import os, sys, json
import numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from T1_graded_blankets import n, net
from T4_observer_dependence import rho_from
rng = np.random.default_rng(5050); srng = np.random.default_rng(6161); dh, du, n_ok = [], [], 0
for g in range(50):
    A = net(rng, 0.0); rng.choice(n, 2, replace=False); Ao = A[:n, :n]; W = Ao * (1 - np.eye(n))
    s = 10 ** srng.uniform(0, 1, n); Sc = lambda X: s[:, None] * X / s[None, :]
    r0n, r3n = rho_from(W ** 2, range(n)), rho_from(expm(3 * Ao) ** 2, range(n))
    r0s, r3s = rho_from(Sc(W) ** 2, range(n)), rho_from(Sc(expm(3 * Ao)) ** 2, range(n))
    am = lambda r: int(np.nanargmax(r))
    dh.append(am(r0s) != am(r3s)); du.append(am(r0n) != am(r0s))
out = dict(fast_vs_slow_differs_random_units=float(np.mean(dh)), fast_noise_vs_random_units_differs=float(np.mean(du)), n=len(dh))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M11.json'), 'w'), indent=1)
