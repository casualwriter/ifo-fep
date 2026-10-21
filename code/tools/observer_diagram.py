"""observer_diagram.py - Figure 1 of 'Which Observer?': the same partition seen by many observers.
Uses blanket_report.py demo systems (Case 1: split blanket, M = 0; Case 2: active state hears outside, M != 0).
Heat map of the leak for observers who read the blanket every delta, over time scale h and read interval delta,
plus the continuous-read limit (delta -> 0, exact Riccati) and the total-effect observer (does not read B)."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from blanket_report import read_leak, cont_leak, total_leak
rng = np.random.default_rng(1)
I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + Ac
hears = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}
A = -1.5 * np.eye(8)
for i in range(8):
    for j in hears[i]:
        if i != j and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
A[4, 0] = 0.5
m = np.linalg.eigvals(A).real.max()
if m > -0.1: A -= (m + 0.1) * np.eye(8)
A2 = A.copy(); A2[4, 6] = 0.3; D = np.eye(8)
hs = [0.25, 0.5, 1.0, 2.0, 4.0]; ks = [1, 2, 4, 8, 16]          # delta = h / k
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2)); fig.subplots_adjust(wspace=0.45)
for ax, (Am, title) in zip(axes, [(A, "Case 1: split blanket (M = 0)"), (A2, "Case 2: active state also hears outside (M ≠ 0)")]):
    Z = np.array([[read_leak(Am, D, I, B, E, h, h / k) for h in hs] for k in ks] + [[cont_leak(Am, D, I, B, E, h) for h in hs]])
    T = [total_leak(Am, I, B, E, h) for h in hs]
    im = ax.imshow(np.log10(np.maximum(Z, 1e-12)), origin="lower", aspect="auto", cmap="viridis", vmin=-12, vmax=0)
    ax.set_xticks(range(len(hs))); ax.set_xticklabels([str(h) for h in hs]); ax.set_xlabel("time scale h")
    ax.set_yticks(range(len(ks) + 1)); ax.set_yticklabels([f"δ=h/{k}" for k in ks] + ["continuous"])
    ax.set_title(title, fontsize=10)
    for j, t in enumerate(T): ax.text(j, len(ks) + 0.6, f"{t:.2f}", ha="center", fontsize=7, color="black")
    ax.set_ylim(-0.5, len(ks) + 1.0)
    ax.text(-0.45, len(ks) + 0.6, "no read:", fontsize=7, ha="right")
    for (r, c), v in np.ndenumerate(Z): ax.text(c, r, "0" if v < 1e-12 else f"{v:.0e}", ha="center", va="center", fontsize=6, color="white")
fig.colorbar(im, ax=axes, label="log10 leak (reads blanket)")
fig.suptitle("Same system, same partition, many observers. Wiring (sparse coupling) says 'exact' in both cases.", fontsize=10)
out = os.path.join(HERE, "..", "..", "paper", "figs", "fig1_observer_diagram.png")
fig.savefig(out, dpi=150, bbox_inches="tight"); print("saved", os.path.abspath(out))
