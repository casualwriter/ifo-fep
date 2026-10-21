"""M10c: start effect for CURRENT states vs PATHS (reconciles Virgo et al. 2022 with M10). Prospectively specified.
Same model class and seed as M10 (split boundary, stationary start not known, spectral radius <= 0.5, 100 models).
Measure A (current states, Virgo's statement): I(s_T; z_T | b_0..b_T). Measure B (paths, Friston 2023 / DMBD statement):
I(s path; z path | b path) (= M10).
Predictions: A fades with history length, median A(T=32)/A(T=8) < 0.2 (ergodic forgetting of the start);
B stays constant, ratio in [0.8, 1.25]."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "M10_dmbd_transient.py"), encoding="utf-8").read().split("\nout = {}")[0]
src = src.replace("    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp))",
                  "    Sc, Zc = [T * n + k for k in s], [T * n + k for k in z]\n"
                  "    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp)), 0.5 * (L(Sc + Bp) + L(Zc + Bp) - L(Bp) - L(Sc + Zc + Bp))")
exec(src)
v = np.array([[cmi(A, True, T) for T in (8, 16, 32)] for A in [gen(True) for _ in range(100)]])   # models x T x (path, current)
out = dict(path_median_T8_16_32=np.median(v[:, :, 0], 0).tolist(), current_median_T8_16_32=np.median(v[:, :, 1], 0).tolist(),
           path_ratio32_8=float(np.median(v[:, 2, 0] / v[:, 0, 0])), current_ratio32_8=float(np.median(v[:, 2, 1] / np.maximum(v[:, 0, 1], 1e-300))))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M10c.json'), 'w'), indent=1)
