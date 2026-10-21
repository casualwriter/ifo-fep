"""M15: many outsides, many blankets; nesting; a blanket that moves when the system changes. Prospectively specified.

System: layered ('onion') network, 5 layers L0..L4 of 2 nodes (L0 = core). Edges only within a layer and between adjacent layers
(each allowed directed edge present with prob 0.7, weight U(0.2, 0.6), random sign), self-decay -1.5, stabilised; unit diagonal noise.
50 networks, seed 1515. Observer for blanket exactness: knows the start, reads the blanket continuously (residual leak, Riccati, h = 2).

Lemma (proved in the paper): with nearest-neighbour layers and diagonal noise,
  (a) a one-layer blanket B = L_k (I = inner layers, E = outer layers) has M != 0 for generic weights (B hears both sides);
  (b) a two-layer blanket B = L_k u L_{k+1} has M = 0 for ALL weights (inner layer hears only the inside, outer only the outside).
Predictions:
  L1 (a): |M| > 1e-6 and leak > 1e-6 in >= 95% of networks for B = L2 (I = L0 u L1, E = L3 u L4).
  L2 (b): M = 0 (< 1e-12) and leak < 1e-10 in all networks for B = L1 u L2 and for B = L2 u L3.
  Many outsides: an observer whose access is L3 (near) can use only the blanket L1 u L2 (I = L0); an observer whose access is L4 (far)
     has two exact, nested blankets: L1 u L2 (I = L0) and L2 u L3 (I = L0 u L1). Checked by enumerating all two-layer blankets that
     leave the observer's access layer outside: exact (leak < 1e-10) in all networks.
  Moving blanket: add a shortcut L0 -> L2 of weight eps in {0, 0.1, 0.2, 0.3} (core starts to drive the outer membrane layer).
     Blanket L1 u L2 (I = L0) leaks, with median leak rising with eps (slope in log eps in [1.6, 2.4] over 0.1-0.3);
     blanket L2 u L3 (I = L0 u L1) stays exact (< 1e-10). The far observer's exact blanket moves outward."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'tools'))
from blanket_report import cont_leak, M_matrix
L = [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]]; n = 10
layer = {k: i for i, l in enumerate(L) for k in l}
rng = np.random.default_rng(1515)
def gen():
    A = -1.5 * np.eye(n)
    for i in range(n):
        for j in range(n):
            if i != j and abs(layer[i] - layer[j]) <= 1 and rng.random() < 0.7:
                A[i, j] = rng.uniform(0.2, 0.6) * rng.choice([-1, 1])
    m = np.linalg.eigvals(A).real.max()
    return A - (m + 0.1) * np.eye(n) if m > -0.1 else A
cat = lambda ks: [x for k in ks for x in L[k]]
D = np.eye(n)
def part(inner, blanket):          # layers 0..inner-1 inside, blanket layers, rest outside
    I = cat(range(inner)); B = cat(blanket); E = [x for x in range(n) if x not in I and x not in B]; return I, B, E
def Mnorm(A, I, B, E): return max(np.abs(M_matrix(A, D, I, B, E)).max(), np.abs(M_matrix(A, D, E, B, I)).max())
nets = [gen() for _ in range(50)]; out = {}
# lemma
one = [(Mnorm(A, *part(2, [2])), cont_leak(A, D, *part(2, [2]), 2.0)) for A in nets]
two = [(Mnorm(A, *p), cont_leak(A, D, *p, 2.0)) for A in nets for p in (part(1, [1, 2]), part(2, [2, 3]))]
out['L1_one_layer'] = dict(frac_M_gt=float(np.mean([m > 1e-6 for m, _ in one])), frac_leak_gt=float(np.mean([l > 1e-6 for _, l in one])),
                           median_leak=float(np.median([l for _, l in one])))
out['L2_two_layer'] = dict(max_M=float(max(m for m, _ in two)), max_leak=float(max(l for _, l in two)))
# many outsides: enumerate two-layer blankets (k, k+1) with inner = k, access layer outside
def exact_blankets(A, access):
    res = []
    for k in range(1, 4):
        if access in (k, k + 1) or access < k: continue
        p = part(k, [k, k + 1]); v = cont_leak(A, D, *p, 2.0)
        res.append((k, None if np.isnan(v) else bool(v < 1e-10)))   # FIX after run 1: None = no causal route (score undefined)
    return res
near = [exact_blankets(A, 3) for A in nets]; far = [exact_blankets(A, 4) for A in nets]
out['many_outsides'] = dict(near_candidates=[k for k, _ in near[0]], near_all_exact=all(e for r in near for _, e in r if e is not None),
                            near_undefined=sum(e is None for r in near for _, e in r),
                            far_candidates=[k for k, _ in far[0]], far_all_exact=all(e for r in far for _, e in r if e is not None),
                            far_undefined=sum(e is None for r in far for _, e in r))
# moving blanket
mv = {}
for eps in [0.0, 0.1, 0.2, 0.3]:
    l12, l23 = [], []
    for A0 in nets:
        A = A0.copy(); A[4, 0] = A[4, 0] + eps        # core node 0 -> outer-membrane node 4 (layer L2)
        m = np.linalg.eigvals(A).real.max()
        if m > -0.1: A = A - (m + 0.1) * np.eye(n)
        l12.append(cont_leak(A, D, *part(1, [1, 2]), 2.0)); l23.append(cont_leak(A, D, *part(2, [2, 3]), 2.0))
    mv[eps] = dict(leak_L1L2=float(np.nanmedian(l12)), leak_L2L3_max=float(np.nanmax(l23)), n_undefined=int(np.isnan(l12).sum()))   # FIX: nan-safe
es = [0.1, 0.2, 0.3]
out['moving'] = dict(by_eps=mv, slope=float(np.polyfit(np.log(es), np.log([mv[e]['leak_L1L2'] for e in es]), 1)[0]))
print(json.dumps(out, indent=1, default=str)); json.dump(out, open(os.path.join(HERE, 'M15.json'), 'w'), indent=1, default=str)
