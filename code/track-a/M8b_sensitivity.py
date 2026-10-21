"""M8b: how the floor's size depends on the observer (referee request C6). Riccati floor (M8 code), Loop 3 nets (seed 9191), eps = 0.1.
Prediction (fixed before running): the median floor changes by more than 2x across h in {0.5, 1, 2, 4}, and by more than 2x across
D_E scale in {0.25, 1, 4}; across push size sigma_u in {0.01, 0.1, 1} it changes by less than 2x (small-push regime up to 1)."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M8_riccati_floor.py"), encoding="utf-8").read().split("\nnets = [gen()")[0])
nets = [gen() for _ in range(200)]; As = []
for M0 in nets: M = M0.copy(); M[4, 6] = 0.1; As.append(stab(M))
def run(h=2.0, s=0.1, de=1.0):
    global H_, su
    H_, su = h, s; D = np.eye(8); D[6, 6] = D[7, 7] = de
    return float(np.median([riccati_floor(A, D) for A in As]))
out = dict(h={h: run(h=h) for h in [0.5, 1, 2, 4]}, sigma_u={s: run(s=s) for s in [0.01, 0.1, 1.0]}, D_E={d: run(de=d) for d in [0.25, 1, 4]})
out['ratios'] = {k: max(v.values()) / min(v.values()) for k, v in out.items()}
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M8b.json'), 'w'), indent=1)
