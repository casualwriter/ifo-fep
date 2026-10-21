"""F5: small pushes. How does the push + sham verification test (F3) scale with push size u and number of pushes K?
Real systems often allow only small pushes (u < 2 noise units).

Fixed before running:
 - networks: F3_verify20.net (20 nodes, I/B/E = 6/6/8), 50 per drive level c in {0, 0.4}, network seed 2020;
   the same networks and observational data for every (u, K) cell.
 - test: F3_verify20.intervention_test (initial slope, Bonferroni 0.05), u in {0.5, 1, 2, 3}, K in {1000, 4000}.
 - outputs: false rejection of the true blanket, detection of the corrupted blanket.

Prediction: the slope signal grows with u and its noise does not, so z ~ u * sqrt(K). Detection depends on u^2 K only:
 (u=1, K=4000) ~ (u=2, K=1000)  and  (u=0.5, K=4000) ~ (u=1, K=1000)   [u^2 K = 4000 and 1000].
 False rejection stays near the F3 level (3-9%) for every u.
Death clauses: (1) scaling fails if, in a matched pair, the two detection rates have non-overlapping 95% Wilson intervals.
 (2) calibration fails at small pushes if false rejection > 15% in any cell."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'causal-partition-v0_6'))
import F3_verify20 as F
from stats import wilson
U, KS, G = [0.5, 1.0, 2.0, 3.0], [1000, 4000], 50
if __name__ == "__main__":
    out = {}; t0 = time.time()
    for c in [0.0, 0.4]:
        F.rng = np.random.default_rng(2020)
        nets = [F.net(c) for _ in range(G)]
        data = [F.simulate(A, 10000, 5) for A, *_ in nets]
        for u in U:
            for K in KS:
                F.rng = np.random.default_rng(int(1000 * u) + K)
                fr, det = [], []
                for (A, true, bad, isbad), X in zip(nets, data):
                    fr.append(F.intervention_test(A, X, true, K=K, u=u))
                    if isbad: det.append(F.intervention_test(A, X, bad, K=K, u=u))
                key = f'c={c}_u={u}_K={K}'
                out[key] = dict(false_reject=wilson(int(sum(fr)), len(fr)), detect=wilson(int(sum(det)), len(det)), n_bad=len(det))
                print(key, out[key], round(time.time() - t0), 's', flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'F5.json'), 'w'), indent=1,
              default=lambda o: o.item() if hasattr(o, 'item') else str(o))
