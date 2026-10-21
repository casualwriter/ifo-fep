"""M13d (POST HOC, after Trace's v0.12 review): the M13c slope (-0.055) is slower than the proven asymptotic rate 2(g_I + g_E) = -0.08,
which a decay BOUND should not allow if the window were asymptotic. Check: same construction as M13c (hidden slow nodes 0 and 6,
self-rate -0.02, both sides), MI at L = 40, 80, 120, 160; natural-log slope over successive windows.
Expectation (no threshold): the slope steepens toward -0.08 at late L if the M13c slope was pre-asymptotic. Seed 1314 stream as M13c."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "M13c_two_sided.py"), encoding="utf-8").read().split("Ls = (0, 5, 10, 20, 40)")[0]
exec(src)
Ls = (40, 80, 120, 160); rows = []
while len(rows) < 30:
    A = gen(); A1 = hide(A.copy(), 0, [1] + S + Ac); A2 = hide(A1.copy(), 6, [7] + S + Ac)
    if stable(A2, 0.01): rows.append(mi_curve(A2, Ls))
v = np.array(rows)
out = {'L': list(Ls), 'median_MI': np.median(v, 0).tolist()}
for a, b in [(0, 1), (1, 2), (2, 3)]:
    ok = (v[:, a] > 1e-300) & (v[:, b] > 1e-300)
    out['slope_%d_%d' % (Ls[a], Ls[b])] = float(np.median((np.log(v[ok, b]) - np.log(v[ok, a])) / (Ls[b] - Ls[a]))) if ok.any() else None
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M13d.json'), 'w'), indent=1)
