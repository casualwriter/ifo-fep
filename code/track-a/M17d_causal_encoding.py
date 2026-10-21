"""M17d (POST HOC, exploratory; after the round-4 review of fep-blanket v0.8). No prediction.
Nova noted that the M17c (b) examples show snapshot encoding without causal screening, not causal encoding.
Question: in the same 30 small-blanket collider networks (M17c b), does the CAUSAL reader (history) have a synchronisation map?
Also counts, over all networks of M17b and M17c, how often the causal reader has a map (R2 > 0.999).
Same generators and seeds as M17b (1718) and M17c (1719)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tools'))
import blanket_report as br
res = {}
# M17c stream: replay part (a) draws, then part (b)
src = open(os.path.join(HERE, 'M17c_sensory_only.py'), encoding='utf-8').read().split("\nout = {}")[0]
g = {'__file__': os.path.join(HERE, 'M17c_sensory_only.py')}; exec(src, g)
I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]
split = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}
c = 0
while c < 30:
    A = g['gen'](8, split)
    if g['stable'](A): c += 1
I6, B6, E6 = [0, 1], [2, 3], [4, 5]
h6 = {0: [0, 1, 2, 3], 1: [0, 1, 2, 3], 3: [0, 1, 2, 3], 2: [2, 3, 4, 5], 4: [2, 3, 4, 5], 5: [2, 3, 4, 5]}
rs, rh = [], []
while len(rs) < 30:
    A = g['gen'](6, h6, extra=[(2, 0), (2, 4)])
    if g['stable'](A):
        rs.append(br.sync_r2(A, np.eye(6), I6, B6, E6, 'snapshot')); rh.append(br.sync_r2(A, np.eye(6), I6, B6, E6, 'history'))
rs, rh = np.array(rs), np.array(rh)
res['small_collider'] = dict(n=30, snapshot_maps=int((rs > 1 - 1e-9).sum()), history_maps=int((rh > 0.999).sum()),
                             history_R2_median=float(np.median(rh)), history_R2_max=float(rh.max()))
m17b = json.load(open(os.path.join(HERE, 'M17b.json')))
res['M17b_history_maps'] = {k: int(round(30 * (1 - m17b[k]['H_frac_lt_0999']))) for k in ['split', 'collider']}
res['causal_maps_total'] = res['small_collider']['history_maps'] + sum(res['M17b_history_maps'].values())
res['networks_total'] = 90
print(json.dumps(res, indent=1)); json.dump(res, open(os.path.join(HERE, 'M17d.json'), 'w'), indent=1)
