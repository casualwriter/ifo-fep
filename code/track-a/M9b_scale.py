"""M9b (post hoc, referee request): scale for M9. Same models (seed 9595) and cases; report the conditional MI
I(s path; z path | b path) relative to the unconditional I(s path; z path), and per time step."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "M9_dmbd_path_factorisation.py"), encoding="utf-8").read().split("\nout = {}")[0]
src = src.replace("    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp))",
                  "    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp)), 0.5 * (L(Sp) + L(Zp) - L(Sp + Zp))")
exec(src)
out = {}
for key, split, stat in [('D1_general_x0known', False, False), ('D2_split_x0known', True, False), ('D3_split_stationary', True, True)]:
    r = np.array([cmi(gen(split), stat) for _ in range(200)])
    out[key] = dict(cond_median=float(np.median(r[:, 0])), uncond_median=float(np.median(r[:, 1])),
                    ratio_median=float(np.median(r[:, 0] / np.maximum(r[:, 1], 1e-300))), cond_per_step_median=float(np.median(r[:, 0]) / 8))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M9b.json'), 'w'), indent=1)
