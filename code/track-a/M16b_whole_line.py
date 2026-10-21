"""M16b (POST HOC, after the prior-art check of 20 Oct 2026; no prediction).
Fasen-Hartmann & Schenk (2025): for a stationary OU process, whole-line path screening iff the inverse-spectral I-E block
K + i w (C' - G) vanishes for all w. Example with no wiring cut but a cancellation (G = C', K = 0): the whole-line block is zero,
while the finite-window path MI on [0, 1] with a known start (Euler dt = 0.02, exact Gaussian) is positive (endpoint term, Prop. 1)."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
A = np.array([[-1, 0, .5], [1, -2, 1], [.5, 0, -1.]]); D = np.eye(3)
blocks = {}
for w in [0, .3, 1, 2.5, 10]:
    g = (-1j * w * np.eye(3) - A).T @ np.linalg.inv(D) @ (1j * w * np.eye(3) - A)
    blocks[str(w)] = abs(g[0, 2])
h, dt = 1.0, 0.02; n = int(round(h / dt)); F = np.eye(3) + dt * A
P = [np.linalg.matrix_power(F, k) for k in range(n + 1)]
C = np.zeros((3 * n, 3 * n))
for a in range(n):
    for b in range(n):
        C[3*a:3*a+3, 3*b:3*b+3] = sum(P[a - k] @ (dt * D) @ P[b - k].T for k in range(min(a, b) + 1))
idx = lambda c: [3 * t + c for t in range(n)]
I, B, E = idx(0), idx(1), idx(2); ld = lambda X: np.linalg.slogdet(X)[1]; S = lambda s: C[np.ix_(s, s)]
mi = 0.5 * (ld(S(I + B)) + ld(S(E + B)) - ld(S(B)) - ld(S(I + B + E)))
out = dict(max_abs_spectral_IE=max(blocks.values()), finite_window_path_MI=float(mi))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M16b.json'), 'w'), indent=1)
