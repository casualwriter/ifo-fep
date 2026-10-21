"""M8d (POST HOC, after M8c missed its 1% tolerance at t = 0.01 with 2.0%): if the leading term -(t^3/3) D_I M D_E is right, the
relative error should shrink about linearly in t. Check at t = 0.01, 0.001, 0.0001 on the same 200 networks."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M8c_nats_necessity.py"), encoding="utf-8").read().split("\nnets = [gen()")[0])
nets = [gen() for _ in range(200)]; out = {}
for t in [0.01, 0.001, 0.0001]:
    e = []
    for M0 in nets:
        M = M0.copy(); M[4, 6] = 0.1; got, pred = nopush_cross(stab(M), t)
        if np.abs(pred).max() > 0: e.append(float(np.linalg.norm(got - pred) / np.linalg.norm(pred)))
    out[str(t)] = float(np.median(e))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M8d.json'), 'w'), indent=1)
