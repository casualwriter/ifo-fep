"""M3b (POST-HOC, after M3): L2a failed -- lc of a wiring cut plateaus at ~0.12 as the B read grain gets finer (m = 16).
Diagnosis to test: Thm 8 is wrong. With the start state and the B path known, p(I path, E path | B path) ~ p(I|B) p(E|B) p(B|I,E):
B's own dynamics depend on both I and E, so conditioning on the B path couples I and E (collider). Prediction from this diagnosis:
 (i) remove all E -> B edges (B hears only I and B): direction I -> E lc -> 0 (median at m = 16 < 0.01);
 (ii) keep E -> B but remove I -> E? (already 0) -- control: general nets, direction I -> E only, plateau stays ~0.1.
Same networks as L2a (seed 21), h = 2."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'causal-partition-v0_6'))
import expC
exec(open(os.path.join(HERE, "M3_loop2.py"), encoding="utf-8").read().split("\nout = {}")[0])   # load Obs only (M3 has no main guard)
expC.rng = np.random.default_rng(21); res = {'general': {m: [] for m in [1, 4, 16, 32]}, 'no_E_to_B': {m: [] for m in [1, 4, 16, 32]}}
for g in range(200):
    A, I, B, E = expC.random_planted(0.0); A6 = A[:6, :6].copy()
    A_ff = A6.copy(); A_ff[np.ix_(B, E)] = 0
    m_ = np.linalg.eigvals(A_ff).real.max()
    if m_ > -0.1: A_ff -= (m_ + 0.1) * np.eye(6)
    for key, M in [('general', A6), ('no_E_to_B', A_ff)]:
        for m in res[key]:
            res[key][m].append(Obs(M, 2.0, m).lc(np.eye(6)[:, I], I, B, E))
out = {k: {m: float(np.nanmedian(v)) for m, v in d.items()} for k, d in res.items()}
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M3b.json'), 'w'), indent=1)
