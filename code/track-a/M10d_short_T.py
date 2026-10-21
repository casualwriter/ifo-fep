"""M10d (POST HOC, after M10c: current-state MI was already ~0 at T = 8, so its ratio was not evaluable).
Same models (seed 9696 stream as M10c), current-state and path MI at T = 1, 2, 4, 8."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M10c_current_vs_path.py"), encoding="utf-8").read().split("\nv = np.array")[0])
v = np.array([[cmi(A, True, T) for T in (1, 2, 4, 8)] for A in [gen(True) for _ in range(100)]])
out = dict(T=[1, 2, 4, 8], path_median=np.median(v[:, :, 0], 0).tolist(), current_median=np.median(v[:, :, 1], 0).tolist(),
           current_q90=np.quantile(v[:, :, 1], 0.9, axis=0).tolist())
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M10d.json'), 'w'), indent=1)
