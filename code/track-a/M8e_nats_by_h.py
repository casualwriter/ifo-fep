"""M8e (referee request): numerator and denominator in nats for each h of M8b (eps = 0.1, Loop 3 nets), so the 24-fold change of the
ratio can be split into numerator and denominator changes. Descriptive (no prediction)."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M8c_nats_necessity.py"), encoding="utf-8").read().split("\nnets = [gen()")[0])
nets = [gen() for _ in range(200)]; As = []
for M0 in nets: M = M0.copy(); M[4, 6] = 0.1; As.append(stab(M))
out = {}
for h in [0.5, 1, 2, 4]:
    H_ = h; v = np.array([[nats(A, I, E), nats(A, E, I)] for A in As])     # nets x dir x (num, den)
    out[str(h)] = dict(num_median=float(np.median(v[:, :, 0])), den_median=float(np.median(v[:, :, 1])),
                       ratio_median=float(np.nanmedian(np.max(v[:, :, 0] / np.where(v[:, :, 1] > 1e-12, v[:, :, 1], np.nan), axis=1))))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M8e.json'), 'w'), indent=1)
