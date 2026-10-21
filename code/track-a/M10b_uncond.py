"""M10b (referee request): conditional AND unconditional I(s path; z path) for each T, same models as M10 (seed 9696). Descriptive."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "M10_dmbd_transient.py"), encoding="utf-8").read().split("\nout = {}")[0]
src = src.replace("    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp))",
                  "    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp)), 0.5 * (L(Sp) + L(Zp) - L(Sp + Zp))")
exec(src); out = {}
for key, split, stat in [('G_general_known', False, False), ('S_split_stationary', True, True), ('K_split_known', True, False)]:
    v = np.array([[cmi(A, stat, T) for T in (8, 16, 32)] for A in [gen(split) for _ in range(100)]])   # models x T x (cond, uncond)
    out[key] = dict(cond_median_T8_16_32=np.median(v[:, :, 0], 0).tolist(), uncond_median_T8_16_32=np.median(v[:, :, 1], 0).tolist())
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M10b.json'), 'w'), indent=1)
