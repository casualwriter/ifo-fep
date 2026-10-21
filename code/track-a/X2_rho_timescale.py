"""X2 (EXPLORATORY, not pre-registered): crossing share seen by an observer with time resolution h (exact, linear).
An observer pushes I (or E) and reads the response after time h. Response matrix R(h) = expm(A h) (observed block).
rho(h) = min over directions of  sum R_BI^2 / (sum R_BI^2 + sum R_EI^2).  h -> 0 gives the coupling-based rho of T1 (1 for a valid
blanket); larger h lets influence pass through B and reach E, so a valid blanket looks leaky to a slow observer.
Networks: the 6-node planted networks of Section 3.4 (expC.random_planted, valid blanket, no hidden drive), 200 networks.
Also: blanket time tau_B = smallest h with rho(h) < 0.5, in units of the nodes' own decay time (1)."""
import os, sys, json
import numpy as np
from scipy.linalg import expm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'causal-partition-v0_6'))
import expC
from expC import random_planted
HS = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
def rho_h(A, I, B, E, h):
    R = expm(A * h)[:6, :6] ** 2
    def one(src, tgt):
        s = lambda t: R[np.ix_(t, src)].sum()
        d = s(B) + s(tgt)
        return s(B) / d if d > 1e-15 else np.nan
    v = [one(I, E), one(E, I)]
    return np.nan if np.isnan(v).all() else float(np.nanmin(v))   # a side with no outgoing influence is skipped
expC.rng = np.random.default_rng(21)
rows = []
for g in range(200):
    A, I, B, E = random_planted(0.0)
    r = [rho_h(A, I, B, E, h) for h in HS]
    fine = np.linspace(0.01, 10, 1000); rf = np.array([rho_h(A, I, B, E, h) for h in fine[::10]])
    if np.isnan(rf).all(): rows.append(dict(rho=r, tau=None)); continue
    tau = float(fine[::10][np.argmax(rf < 0.5)]) if (rf < 0.5).any() else None
    rows.append(dict(rho=r, tau=tau))
R = np.array([x['rho'] for x in rows]); R = R[~np.isnan(R).any(1)]; taus = [x['tau'] for x in rows]
out = {'h': HS, 'networks_used': len(R), 'rho_median': np.median(R, 0).round(3).tolist(), 'rho_q10': np.quantile(R, 0.1, 0).round(3).tolist(),
       'rho_q90': np.quantile(R, 0.9, 0).round(3).tolist(),
       'tau_median': float(np.median([t for t in taus if t is not None])), 'tau_none': sum(t is None for t in taus)}
print(json.dumps(out))
json.dump(out, open(os.path.join(HERE, 'X2.json'), 'w'), indent=1)
