# Finite-Window Path Screening in Linear Diffusions: An Exact Condition and a Residual-Leak Formula

**Authors:** CK Hung, Emmy Team (Amoy Technology)  
**Date:** 22 October 2026  
**Version:** 1.0.1 (October 2026)  
**License:** CC BY 4.0 (code: MIT)  
**DOI:** https://doi.org/10.5281/zenodo.23524008 (all versions; version 1.0 is 10.5281/zenodo.23524009)  
**Main paper:** "Interface With Whom? Markov Blanket Claims Depend on the Observation Protocol" (version 1.0, published together with this note).  

## Abstract

A wiring cut — no direct coupling between internal and external states — does not by itself make internal and external paths conditionally independent given the blanket path. For linear diffusions with noise that is independent across the three groups, and an observer who knows the start state and reads the blanket continuously, we show that screening holds if and only if there is a wiring cut and
$$
M = A_{BI}^{\top} D_B^{-1} A_{BE} = 0 ,
$$
where $A$ is the drift and $D$ the noise covariance.

A split blanket — sensory states driven by the outside, active states driven by the inside, with uncorrelated noise between the two kinds — guarantees $M = 0$. Its sufficiency is stated in the path-integral formulation of the free energy principle. In discrete time, the analogous result for current states given the sensorimotor history follows from d-separation on the time-unrolled loop, and a continuous-time analogue was expected.

We give that continuous-time, path-level result, and add seven things:
- the noise-weighted algebraic condition on a finite window, with a necessity proof that also covers the wiring cut (for a stationary reader of the whole time line, a spectral form of the condition is known, and there the wiring cut is not necessary);
- a long-history result: an observer that has read the blanket long enough forgets the start, because of detectability, at a rate set by the slow modes that both sides show only poorly through the blanket;
- a causal-reader result: a reader of the past blanket only, with no future readings, sees the current internal and external states screened off if and only if $M = 0$; in our networks a reader of the sensory states alone did not;
- a separation of screening from encoding: a synchronisation map from internal to external expectations needs a different condition; by example, neither condition implies the other (for encoding read from the current blanket), and for the causal reader we found no exact map at all;
- an explicit condition on the start state, which matters for paths even when it fades for current states;
- an exact Riccati formula for the residual leak when $M \neq 0$;
- a check that this residual leak vanishes when two different violations cancel in $M$.

Dynamic Markov Blanket Detection justifies path factorisation by the absence of direct interaction alone, and its code lets a boundary latent depend on both sides. For generic parameters in that model class, the factorisation fails, and the leftover dependence grows with path length. With a split boundary, it fails only through the start state.

When screening fails, the measured size depends on the system and on the observation protocol together. For example, the residual leak ratio changes about 24-fold across time scales. In layered systems with diagonal blanket noise, a two-layer shell has $M = 0$ for all weights. In one real-data case (fMRI effective connectivity), the present estimator and calibration do not support a decision.

## 1. Introduction

In the free energy principle (FEP), a particle is defined by a Markov blanket: internal and external states are conditionally independent given blanket states. In the path-integral formulation, this is stated for **paths** (Friston et al., 2023).

Two lines of work examine the condition:
- Sparse coupling does not in general give conditional independence of states at steady state (Aguilera et al., 2022; Heins & Da Costa, 2022).
- Conditional independence given the current blanket state is rarely expected. In discrete time it holds given the blanket's history and the initial states, by d-separation on the time-unrolled sensorimotor loop; an analogous continuous-time result is expected (Virgo, Rosas & Biehl, 2022).

Detection methods now use the path-level condition directly (Beck & Ramstead, 2025).

This paper makes the path-level condition exact for linear diffusions and computes what happens when it fails:
1. **Proposition 1.** A wiring cut is necessary but not sufficient. For an observer who knows the start state, the extra condition is $M = 0$. For an observer who does not, a condition on the start law is added.
   **Proposition 1b.** With a long history of blanket readings, the start condition is not needed: the start is forgotten through the blanket.
   **Proposition 1c.** A reader of the past blanket only (no future readings) needs the same condition, $M = 0$, for the current states. In our networks it also needed the active states: a reader of the sensory states alone saw a leak even when $M = 0$.
   **Section 3.2d.** Screening is not encoding: the FEP's synchronisation map needs a separate condition.
2. **Proposition 2.** When $M \neq 0$, the residual leak of a push observer who reads the blanket continuously is given exactly by a Riccati equation. Its size depends on the remaining dynamics, the noise and the protocol, not on $M$ alone.
3. **Section 4.** The stated justification of the path condition in DMBD does not mention $M = 0$ and the start condition, which our results show are needed. We say this with respect for that work and its code, which we use.

Section 5 then maps how the observation protocol changes what is seen: time-scale laws, a slow-observer limit, tests of whether the best partition depends on the time scale, and a slow-mode control.

**Scope.** Linear Gaussian dynamics (the local linearisation of smooth systems near a fixed point), small networks, and exact computation. Nonlinear systems are open. One real-data case (Section 5.7) is inconclusive.

**Evidence labels.**
- **Algebra check.** A numerical confirmation of a derived statement. Several of our checks are of this kind; they verify the theorem and the code, but they are not independent empirical evidence.
- **Prospectively specified.** The prediction and the failure condition were written in the script header before the run. There is no separate time-stamped record, so we do not count these as pre-registered confirmations.
- **Post hoc.** The analysis was chosen after seeing results.

## 2. Set-up

### 2.1 System
$$
dx = A x\, dt + D^{1/2} dW . \tag{1}
$$
- $A$ is stable. The states are split into internal $I$, blanket $B$ and external $E$.
- **Wiring cut:** $A_{IE} = A_{EI} = 0$.
- **Noise:** $D$ is block-diagonal over $(I, B, E)$, with each block positive definite.
- **Noise units:** a representation in which $D = I$. All share-type scores below are computed in noise units, unless a variation of $D$ is stated.
- **Split blanket:** $B = S \cup A_c$, sensory and active, such that:
  - sensory rows do not depend on $I$ ($A_{SI} = 0$);
  - active rows do not depend on $E$ ($A_{A_c E} = 0$);
  - $D_B$ is block-diagonal over $(S, A_c)$.

### 2.2 Three groups of parameters

| Group | Parameters |
|---|---|
| System | $A$, $D$, and the physical distribution of start states |
| Observation protocol | time scale $h$; read interval $\delta$ on the blanket; push and read sets; push distribution; knowledge of the start state |
| Representation | units and normalisation of a score |

Proposition 1 holds within one protocol class: the start state is known and the blanket is read continuously. Within that class, whether screening is exact is a property of the system ($M$). Outside it, the protocol also matters. How large a failure looks depends on the system and the protocol together.

### 2.3 Pushes and scores
**Push.** A push acts on all internal states at time 0:
$$
x(0^+) = x(0^-) + \sigma_u P_I u, \qquad u \sim N(0, I_{|I|}). \tag{2}
$$
- $P_I$ is the injection of the internal coordinates into the full state.
- "Start state known" means that $x(0^-)$ is known.

**Total-effect share.** This score is for an observer who does not read the blanket between push and read. With $R(h) = e^{Ah}$:
$$
\ell_{I\to E}(h) = \frac{\sum_{e,i} R_{ei}^2}{\sum_{b,i} R_{bi}^2 + \sum_{e,i} R_{ei}^2}, \qquad \rho(h) = 1 - \max(\ell_{I\to E}, \ell_{E\to I}). \tag{3}
$$

**Conditional push-information ratio (CPIR).** This score is for an observer who reads the blanket every $\delta$:
$$
\ell^{c}_{\delta}(h) = \frac{I\big(u; x_E(h) \mid x(0^-), x_B(\delta), \dots, x_B(h)\big)}{I\big(u; x_E(h) \mid x(0^-)\big)} . \tag{4}
$$
The **residual leak** is the continuous-read limit $\ell^c_0(h) = \lim_{\delta \to 0} \ell^c_\delta(h)$. We use $\rho^c = 1 - \max$ over the two directions.

Properties of the CPIR:
- **Units.** It is invariant to rescaling when the push is transformed with the variables (the same experiment), because it is built from mutual information (Section 2.4).
- **Not a causal share.** The numerator includes explaining away, so the ratio can exceed 1.
- **Domain.** It is defined only when the denominator is positive, that is, when the push has a causal route to $E$.
  - Counter-example: a blanket that is a sink driven by both sides has a zero denominator but can have a positive numerator.
  - In our 8-node networks, 9 of 400 directions had a denominator below $10^{-12}$; they were excluded.
- **Not a proven lower bound.** We do not claim that $\ell^c_0$ bounds $\ell^c_\delta$ from below. More readings can add collider information. In our tests $\ell^c_\delta$ decreased as $\delta$ fell.
- **Reporting.** We report the numerator in nats alongside the ratio.

### 2.4 Units
A share built from squared couplings changes under $x \to Sx$ with $S$ diagonal.
- **Test:** in 50 random 12-node networks, a random $S$ (log-uniform over one decade).
- **Result:** the best partition by direct couplings changed in 50% of networks. The CPIR changed by at most $5\times10^{-13}$.
- **What is invariant.** In this test the push was transformed together with the variables, so it is the same physical experiment; mutual information is then invariant. If the push is instead redefined as isotropic in the new units, it is a different observer, and its CPIR changes (in one 8-node network, 4.40 vs 1.58). The push distribution is part of the observer (Section 2.2).

## 3. The exact condition and the residual leak

### 3.1 Proposition 1 (path screening)

**Assumptions:**
- $D$ block-diagonal over $(I,B,E)$, with $D_I, D_B, D_E \succ 0$;
- the start state $x(0^-)$ is known.

**Claim:** on $[0,h]$, the internal path and the external path are conditionally independent given the blanket path **if and only if**
$$
A_{IE} = 0,\qquad A_{EI} = 0 \quad\text{(wiring cut)}, \qquad\text{and}\qquad M \equiv A_{BI}^{\top} D_B^{-1} A_{BE} = 0. \tag{5}
$$

**Corollary.**
- A split blanket (Section 2.1) with a wiring cut gives screening for any weights.
- Suppose $D_B$ is diagonal and the weights are in general position, meaning there is no exact cancellation in the sum $\sum_b d_b^{-1} A_{bI}^{\top}A_{bE}$. Then $M = 0$ requires that no blanket state be driven by both sides.

*Proof.*

**Step 1: the path density.**
- Take the product of the driftless laws of the three groups, each started at its known start point, as the reference law.
- Relative to it, the law of the path on $[0,h]$ has the density $\exp(\ell)$, where
$$
\ell = \sum_{k \in \{I,B,E\}} \Big[ \int_0^h f_k^{\top} D_k^{-1}\, dx_k - \tfrac12 \int_0^h f_k^{\top} D_k^{-1} f_k\, dt \Big],
$$
  and $f_k$ is the drift of group $k$, for example $f_I = A_{II}x_I + A_{IB}x_B + A_{IE}x_E$.
  - Absolute continuity holds for linear drifts by Beneš's condition (Karatzas & Shreve, 1991, Cor. 3.5.16). The Feldman–Hájek theorem for Gaussian measures gives the same result.
  - Novikov's condition, cited in earlier drafts, is not needed and does not hold on every interval.
  - The squared form $-\tfrac12\int|dx_B - f_B\,dt|^2_{D_B^{-1}}$ used informally elsewhere is a formal action. The proof uses only the Itô form above.
- The reference law is a product, so the conditional density of $(x_I, x_E)$ given the blanket path is proportional to $\exp(\ell)$ with $x_B$ held fixed.
- The terms of $\ell$ that are linear in $dx_B$ involve one side only. Since $D_{IB} = D_{EB} = 0$, the covariations vanish. After integration by parts given $x_B$, these integrals are Wiener integrals, defined for almost every blanket path.

**Step 2: the mixed part.**
- **The mixed terms.** The terms of $\ell$ that involve both $x_I$ and $x_E$ are:
$$
\int (A_{IE}x_E)^{\top} D_I^{-1} dx_I + \int (A_{EI}x_I)^{\top} D_E^{-1} dx_E - \int x_I^{\top}\big(A_{II}^{\top}D_I^{-1}A_{IE} + A_{EI}^{\top}D_E^{-1}A_{EE} + M\big)x_E\,dt .
$$
- **The mixed second difference.** Shift the paths along Cameron–Martin directions $\eta_I$ and $\eta_E$. These are absolutely continuous with square-integrable derivative and start at zero, as the known start requires. Integrate the first term by parts. The mixed second difference of $\ell$ is then
$$
\Phi(\eta_I,\eta_E) = \int_0^h \eta_I^{\top}(G - C')\,\dot\eta_E\,dt \;+\; \eta_I(h)^{\top} C'\, \eta_E(h) \;-\; \int_0^h \eta_I^{\top} K\, \eta_E\,dt,
$$
  with $G = A_{EI}^{\top}D_E^{-1}$, $C' = D_I^{-1}A_{IE}$ and $K = A_{II}^{\top}D_I^{-1}A_{IE} + A_{EI}^{\top}D_E^{-1}A_{EE} + M$.

**Step 3: sufficiency.**
- With a wiring cut and $M = 0$, all mixed terms vanish.
- The conditional density is then a product of a function of the internal path and a function of the external path.
- With the start known, the paths are conditionally independent.

**Step 4: necessity.**
- **Φ must vanish everywhere.**
  - Conditional independence requires the conditional density to factorise almost surely, relative to the product reference.
  - Shifts along Cameron–Martin directions give equivalent measures (Cameron–Martin theorem). So the factorisation carries over to shifted paths, and $\Phi$ must vanish for all $\eta_I$ and $\eta_E$.
  - The functional $\ell$ is defined up to null sets. For Cameron–Martin $\eta$, the identity $\ell(x+\eta) = \ell(x) + (\text{linear in }\eta) + (\text{quadratic in }\eta)$ holds almost surely, and the argument uses only this.
- **First, directions with $\eta_I(h) = 0$.**
  - $\int_0^h \eta_I^{\top}\big[(G - C')\dot\eta_E - K\eta_E\big]dt = 0$ for all such $\eta_I$. Hence $(G - C')\dot\eta_E(t) = K\eta_E(t)$ for all $t$.
  - With $\eta_E(t) = b\,t$, this gives $(G - C')b = K b\, t$ for all $t$ and all $b$. So $K = 0$ and $G = C'$.
- **Then, free endpoints.**
  - What remains is $\Phi = \eta_I(h)^{\top}C'\eta_E(h)$, and it must vanish for all endpoint values. So $C' = 0$, that is $A_{IE} = 0$.
  - Then $G = C' = 0$, that is $A_{EI} = 0$.
  - Finally $K = M = 0$.
- **What the endpoint term does.** It rules out the cancellation $A_{EI}^{\top}D_E^{-1} = D_I^{-1}A_{IE}$ raised in review: that cancellation removes the integral term but not the endpoint term.
- **Unknown start.** If the start is not known, and its law has a positive density (for example a non-degenerate Gaussian), the conditional density has the extra factor $p(x_I(0), x_E(0) \mid x_B(0))$. Shifts with $\eta(0) = 0$ do not change it, so it cancels in $\Phi$. The necessity of the wiring cut and of $M = 0$ therefore holds without a known start; Section 3.2 then adds the start condition. $\square$

*Check of the endpoint argument (M16, prospectively specified; Euler discretisation, $dt = 0.02$, $h = 1$, exact Gaussian path MI, 40 networks per case):*

| Case | Prediction | Result (MI in nats) |
|---|---|---|
| Wiring cut | < $10^{-10}$ | $\le 1.4\times10^{-12}$ |
| Symmetric direct couplings $A_{EI} = cR$, $A_{IE} = cR^{\top}$ ($G = C'$) | > $10^{-4}$ in all | min 0.0034, median 0.011 |
| One-way direct coupling | > $10^{-4}$ in all | min 0.0023 |

**A local illustration (filtering, no push).** Condition on the blanket path up to time $t$, with a wiring cut and $x(0^-)$ known. Then $\Sigma_{II}(t) = D_I t + O(t^2)$ and $\Sigma_{EE}(t) = D_E t + O(t^2)$, and the Riccati equation of Proposition 2 gives
$$
\Sigma_{IE}(t) = -\tfrac{t^3}{3} D_I M D_E + O(t^4).
$$
- This is a filtering quantity (blanket read up to $t$), not the path statement itself.
- **Algebra check:** the relative error is 2.0%, 0.20% and 0.020% at $t = 10^{-2}, 10^{-3}, 10^{-4}$ in 200 networks, so it shrinks in proportion to $t$.
- With a push, the leading term of $\Sigma_{uE}$ is of order $t^2$.

### 3.2 Start state
Suppose the start state is not known, and there is a wiring cut and $M = 0$. The conditional law of the paths given the blanket path then contains the factor $p(x_I(0), x_E(0) \mid x_B(0))$. All other factors split by side and are positive.

So, for a start law with a positive density, the paths are conditionally independent if and only if the start law factorises given $x_B(0)$ (all other factors are positive and one-sided, so a product form for almost every blanket path forces this factor to factorise):
- this is necessary and sufficient when $M = 0$;
- at a non-equilibrium steady state it generally fails (Aguilera et al., 2022).

**Current states vs paths.** We checked this in the discrete-time model class of Section 4.2 (post hoc for the short histories): split boundary, stationary start not known, spectral radius ≤ 0.5, 100 models.
- **Current states.** The dependence between the current states given the blanket history fades quickly with history length. Median $I(s_T; z_T \mid b_{0:T})$ = $1.0\times10^{-5}$, $1.4\times10^{-8}$, $8.7\times10^{-14}$ and 0 nats at $T$ = 1, 2, 4, 8.
- **Paths.** The dependence between the paths on $[0, T]$ stays at 0.0048 nats for every $T$. Here the path begins at the unknown start, and the record has no history before the path; a longer $T$ only makes the path longer. With a history *before* the path, the start is forgotten (Proposition 1b).
- **Reading.** For an ergodic system, the start condition matters little for statements about current states with a long history, as Virgo et al. (2022) note. It matters fully for path-level statements, as in Friston et al. (2023) and DMBD.

Friston et al. (2023) write path densities that score fluctuations from a given point, which in effect conditions on the start. The start condition therefore clarifies a convention; it is not a gap in their argument. It matters for methods that place a joint prior on the start state (Section 4).

### 3.2b Proposition 1b (long history: the start is forgotten through the blanket)

**Assumptions:**
- a wiring cut and $M = 0$;
- $D$ block-diagonal with $D_I, D_B, D_E \succ 0$;
- $A$ is stable and the process is stationary;
- the observer does not know the start, but has read the blanket on $[-L, h]$.

**Claim.**
- $I(\text{internal path on }[0,h];\ \text{external path on }[0,h] \mid \text{blanket on }[-L,h]) \to 0$ exponentially as $L \to \infty$, for the stationary (Gaussian) start. In fact it **equals** the mutual information of the start states $x_I(0)$ and $x_E(0)$ given the blanket on $[-L,h]$ (step iv).
- The cross-covariance decays at rate $\gamma_I + \gamma_E$, where $\gamma_I$ and $\gamma_E$ are the decay rates (minus the spectral abscissas) of the stationary closed-loop matrices $C_I = A_{II} - P_I G_{II}$ and $C_E = A_{EE} - P_E G_{EE}$, with $G_{II} = A_{BI}^{\top}D_B^{-1}A_{BI}$ and $G_{EE} = A_{BE}^{\top}D_B^{-1}A_{BE}$. Precisely: $\|\Sigma_{IE}(L)\| \le c_\epsilon\, e^{-(\gamma_I + \gamma_E - \epsilon)L}$ for every $\epsilon > 0$. This is an asymptotic rate; defective closed-loop matrices add polynomial factors. For small dependence the mutual information is about $\tfrac12\sum_k \rho_k^2$ (canonical correlations), so it decays at about twice this rate.
- If $M \neq 0$, the dependence does not vanish; with Proposition 1 this also shows that the wiring cut is necessary in the long-history limit. By the chain rule it is at least $I\big(x_I(0,h];\, x_E(0,h] \mid x(0), x_B[0,h]\big)$. By the Markov property this does not depend on $L$, and it is positive by Proposition 1.

*Proof.*

**(i) Detectability (PBH test).**
- Let $F = \mathrm{diag}(A_{II}, A_{EE})$ and $H = [A_{BI}\ A_{BE}]$.
- Suppose $(F - \lambda)v = 0$ and $Hv = 0$ for some $v = (v_I, v_E) \neq 0$. With the wiring cut,
$$
A\,(v_I, 0, v_E)^{\top} = (A_{II}v_I,\; A_{BI}v_I + A_{BE}v_E,\; A_{EE}v_E)^{\top} = \lambda\,(v_I, 0, v_E)^{\top}.
$$
- So $\lambda$ is an eigenvalue of $A$, and $\mathrm{Re}\,\lambda < 0$. Hence $(F, H)$ is detectable. The cases $v_E = 0$ and $v_I = 0$ show that $(A_{II}, A_{BI})$ and $(A_{EE}, A_{BE})$ are detectable too.
- The interior alone may be unstable. Its unstable modes are then visible through the blanket.

**(ii) Convergence.**
- $Q = \mathrm{diag}(D_I, D_E) \succ 0$, so $(F, Q^{1/2})$ is stabilisable.
- Given the blanket path (the terms $A_{IB}x_B$ and $A_{EB}x_B$ are known inputs), the filtering covariance of $(x_I, x_E)$ solves the equation below. The inputs depend on the observed path itself; the equation still applies, because the joint process is Gaussian and the conditional covariance does not depend on the observed values:
$$
\dot\Sigma = F\Sigma + \Sigma F^{\top} + Q - \Sigma G\Sigma, \qquad G = H^{\top}D_B^{-1}H .
$$
- From any positive semi-definite start, it converges exponentially to the unique stabilising solution $\bar P$ of the algebraic Riccati equation (Anderson & Moore, 1979, Ch. 4).

**(iii) The limit is block-diagonal.**
- With $M = 0$, $G = \mathrm{diag}(G_{II}, G_{EE})$.
- Let $P_I$ and $P_E$ be the stabilising solutions of the two separate equations; they exist by (i).
- Then $\mathrm{diag}(P_I, P_E)$ solves the full equation, and its closed loop $\mathrm{diag}(C_I, C_E)$ is stable. By uniqueness, $\bar P = \mathrm{diag}(P_I, P_E)$.
- **Rate.** The cross block obeys the linear equation $\dot\Sigma_{IE} = (A_{II} - \Sigma_{II}G_{II})\Sigma_{IE} + \Sigma_{IE}(A_{EE} - \Sigma_{EE}G_{EE})^{\top}$, which is homogeneous in $\Sigma_{IE}$. Its coefficients converge to $C_I$ and $C_E$. So $\Sigma_{IE} \to 0$, with the asymptotic rate $\gamma_I + \gamma_E$ stated in the Claim.

**(iv) From covariance to paths.**
- **Factorisation given the start.** Given the blanket on $[-L,h]$, the conditional law of the paths on $[0,h]$ factorises given $(x_I(0), x_E(0))$. This is Proposition 1 on $[0,h]$, with $M = 0$ and a wiring cut.
- **Equal to the start dependence.** Let $W$ be the blanket on $[-L,h]$, $U = (x_I(0), x_I(0,h])$ and $V = (x_E(0), x_E(0,h])$. Given $W$, the conditional density is proportional to $p_f(x_I(0), x_E(0) \mid x_B[-L,0])\, g(U, x_B)\, k(V, x_B)$ (Proposition 1 with $M = 0$, and the Markov property). By the chain rule,
$$
I(U;V \mid W) = I(x_I(0); V \mid W) + \underbrace{I(x_I(0,h]; V \mid x_I(0), W)}_{=0} = I(x_I(0); x_E(0) \mid W) + \underbrace{I(x_I(0); x_E(0,h] \mid x_E(0), W)}_{=0} .
$$
  Each marked term is zero because, given the conditioned variables, the density is a product of a one-sided factor and the rest. So the path MI equals the start MI given the blanket. (Earlier drafts cited the data-processing inequality here; it does not give this bound, because $x_I(0)$ is part of $U$.)
- **Smoothing adds no cross term.** Readings after time 0 multiply the start law only by one-sided factors. So the cross block of the smoothed precision at time 0 equals that of the filtered precision. The filtered precision is $\Sigma(L)^{-1}$, and $\Sigma_{II}$ and $\Sigma_{EE}$ stay away from singular because $P_I, P_E \succ 0$.
- **Conclusion.** The Gaussian MI is continuous in the precision and vanishes when its cross block vanishes. It therefore tends to zero, with a rate set by (iii). $\square$

**Checks** (prospectively specified unless marked):

| Check | Prediction | Result |
|---|---|---|
| M13, discrete time: split blanket, stationary start, history $L$ before the path (100 models) | path MI falls to about 0 | median 0.0048, $1.2\times10^{-8}$, $8\times10^{-14}$, 0, 0 nats at $L$ = 0, 2, 4, 8, 16 |
| M13, blanket state hears both sides ($M \neq 0$) | does not fall | 0.041 → 0.034 (constant from $L$ = 2) |
| M13b (a), continuous time: interior alone unstable, $A$ stable (30 networks) | MI($L$=40) < $10^{-6}$ MI(0) | max ratio $3\times10^{-14}$ |
| M13b (b), random split networks (30) | same | max ratio $5\times10^{-14}$ |
| M13b (c), one interior mode hidden from the blanket, rate −0.02 (30) | forgetting slow (ratio > $10^{-3}$) | **failed**: forgetting stayed fast (MI ≈ $10^{-13}$ at $L$ = 10) |
| M13c (post hoc), hidden slow modes on both sides (30) | slow; natural-log slope of the MI in [−0.16, −0.04] | ratio at $L$ = 40: 0.11; slope −0.055. The proven asymptotic rate is about $-2(\gamma_I+\gamma_E) = -0.08$. The $L \in [10, 40]$ window was pre-asymptotic: M13d (post hoc, same networks, later windows) gives slopes $-0.071$ ($L$ = 40–80), $-0.078$ (80–120) and $-0.080$ (120–160), converging to $-0.08$ |
| M13c (post hoc), hidden slow mode on one side only (30) | fast | ratio about 0 |

**Reading of the failed prediction.** (c) failed because the cross block decays at the **sum** of the two sides' rates (step iii), so one fast side is enough.
- **The bound is in terms of the closed-loop matrices $C_I$ and $C_E$.** Stable modes that the blanket does not show can also forget the start through their own dynamics and noise, so the start is not forgotten *only* through the blanket.
- **Memory of the start survives only if both the inside and the outside have a slow closed-loop mode, that is, a slow mode that the blanket shows only poorly** (strictly hidden modes are the extreme case).
- This is the situation in the example of Virgo et al. (2022): a phone number kept in memory, and also on a note put away in a box. The rate above says how long such shared memory lasts.

### 3.2c Proposition 1c (a causal reader: no future readings needed)

Proposition 1 and 1b condition on the blanket during $[0,h]$. For a reader that is part of the system, that is a smoothing quantity: it uses readings after the time it asks about. This section asks what a **causal** reader needs: one that, at time $t$, has read only the past blanket $(-\infty, t]$.

**Assumptions:** a wiring cut; $D$ block-diagonal with each block positive definite; $A$ stable; stationary process. The reader reads the whole blanket, sensory and active states.

**Claim.** $I\big(x_I(t);\, x_E(t) \mid \text{blanket on } (-\infty, t]\big) = 0$ if and only if $M = 0$.

*Proof.* The stationary filtering covariance of $(x_I(t), x_E(t))$ is the stabilising solution $\bar P$ of
$$
F\bar P + \bar P F^{\top} + Q - \bar P G \bar P = 0 ,
$$
with $F$, $Q$, $G$ as in Proposition 1b; it exists and is unique by detectability (Proposition 1b, step i). The conditional covariance given the infinite past is the limit, as $L \to \infty$, of that given $x_B[-L,t]$ (martingale convergence; for Gaussian processes the conditional covariance does not depend on the observed values), and this limit is $\bar P$. The $I$–$E$ block of $G$ is $M$.
- If $M = 0$, $\bar P$ is block-diagonal (Proposition 1b, step iii).
- If $\bar P_{IE} = 0$, the $I$–$E$ block of the equation reads $-\bar P_I M \bar P_E = 0$. Since $Q \succ 0$, $\bar P_I, \bar P_E \succ 0$, so $M = 0$. $\square$

**What this answers.** A reader of the past blanket only needs the same condition as the path observer. Future readings are not needed for screening of the current states. The path statement on $[0,h]$ itself still conditions on readings during $[0,h]$.

**Scope.** The wiring cut is assumed, not shown necessary, for this reader. Without it, the cross block reads $A_{IE}\bar P_E + \bar P_I A_{EI}^{\top} - \bar P_I M \bar P_E = 0$, which can hold by cancellation.

**The reader needs its own actions.** The interior drives the active states. Does it need to read them? A reader of the sensory states only must also estimate the active states. They carry the interior's influence to the outside ($I \to A \to E$) and are not read, so they open a path around the sensory states. *Check M17c (a)* (prospectively specified; 30 split networks, $M = 0$): the sensory-only reader's dependence $I(x_I(t); x_E(t) \mid$ past sensory states$)$ was above $10^{-4}$ nats in all 30 (median 0.015, range 0.0011–0.23); the whole-blanket reader's was 0 in all. So, in these networks, an exact blanket for the interior needed access to the active states as they are, noise included. This is what d-separation predicts (the path $I \to A \to E$ is open). It is a numerical confirmation, not a statement about physical interiors; special systems may screen with sensory access alone.

*Check M17* (prospectively specified; 8-node networks of Section 3.5, wiring cut, unit noise; 30 split and 30 collider networks; the reader conditions on the whole blanket, sensory and active):

| Measure | Prediction | Result |
|---|---|---|
| Split ($M = 0$): $I(x_I(t); x_E(t) \mid$ past blanket$)$ | < $10^{-10}$ in all | max $2\times10^{-16}$ |
| Collider ($M \neq 0$): same | > $10^{-4}$ in all | min $3.8\times10^{-4}$, median $3.1\times10^{-3}$ nats |
| Identity used in the proof | error < $10^{-9}$ | 0 |

These are algebra checks.

### 3.2d Screening is not encoding

In the FEP, a particle's internal states are also said to *encode* beliefs about external states: the expected internal state given the blanket, $\eta = E[x_I \mid \text{blanket}]$, is said to determine the expected external state, $\mu = E[x_E \mid \text{blanket}]$, through a map $\mu = \sigma(\eta)$, the synchronisation map. Exact screening does not imply this.

**Condition.** In a linear-Gaussian system both expectations are linear in the blanket data, so a map exists if and only if $\mu$ is linear in $\eta$ (a $\sigma(\eta)$-measurable variable that is jointly Gaussian with $\eta$ is almost surely affine in $\eta$). For the snapshot reading (current blanket state only), $\eta = W_I x_B$ and $\mu = W_E x_B$ with $W_k = \Sigma_{kB}\Sigma_{BB}^{-1}$, and, since $\Sigma_{BB} \succ 0$, a map exists if and only if $\ker W_I \subseteq \ker W_E$. If the blanket has more dimensions than the interior, $W_I$ has a kernel, and for generic parameters $W_E$ does not vanish on it. This is the existence condition of Da Costa et al. (2021, Lemma 2.1, with internal and external roles named $\mu$ and $\eta$ the other way round), who also note that it holds generically when the blanket is no larger than the interior. This condition is different from $M = 0$.

**Measure.** $R^2 = 1 - \mathrm{tr}\,\mathrm{Cov}(\mu - E[\mu \mid \eta]) / \mathrm{tr}\,\mathrm{Cov}(\mu)$, the share of the variance of $\mu$ that is linear in $\eta$; a map exists if and only if $R^2 = 1$. We compute it for the snapshot reading and for the causal reader of Proposition 1c (history $(-\infty, t]$, Kalman filter).

*Check M17b* (prospectively specified; same network classes as M17, 30 + 30 networks; positive control: 6-node networks with a 2-state blanket, one sensory and one active, where $W_I$ is generically invertible and so the snapshot map exists; death clause: if the map generically exists, this section is wrong):

| Case | Prediction | Snapshot $R^2$ | History $R^2$ |
|---|---|---|---|
| Control, blanket no larger than the interior | $R^2 = 1$ | min $1 - 10^{-13}$ | — |
| Split ($M = 0$) | $R^2 < 0.999$ in ≥ 90% | 100%; median 0.80 (0.35–0.99) | 100%; median 0.55 (0.12–0.93) |
| Collider ($M \neq 0$) | same | 100%; median 0.73 (0.10–0.97) | 100%; median 0.51 (0.06–0.96) |
| M17c (b): collider ($M \neq 0$, min 0.16), blanket no larger than the interior, 30 | $R^2 = 1$ in all | 29 of 30 (failed as stated); the exception has $\mathrm{rank}\,W_I = 1$, so the kernel condition fails, as the theory says | — |

Exploratory (M17, no threshold): the share of $E[x_E(t) \mid$ past blanket$]$ that is linear in the internal *state* $x_I(t)$, rather than in its expectation, has median 0.22 (split) and 0.19 (collider).

**Reading.**
- **Two conditions, not one.** Screening needs $M = 0$. Encoding through a synchronisation map needs a kernel condition on the blanket's information. Neither implies the other, by example: split networks screen but do not encode (M17b), and the M17c (b) collider networks encode (snapshot reading) but do not screen (causal or path reading; causal-reader dependence at least $3.3\times10^{-4}$ nats). Snapshot screening is a different matter: it can hold with $M \neq 0$ (the stationary precision of Section 5.5 can have a zero $I$–$E$ block while the path reader leaks).
- **No exact causal encoding found.** In the same 30 small-blanket collider networks, the causal reader had no map (median $R^2$ 0.71, max 0.993), and none of the 60 M17b networks had one either: 0 of 90 (check M17d, post hoc). We conjecture that exact causal encoding is non-generic, because the reader's estimate of the two sides generically has full rank. Encoding is then graded, like screening.
- **More history, less encoding** (in these samples). A reader with a long history knows more about both sides, and the two expectations become less tied to each other (median 0.80 → 0.55). We do not claim this as a general law.
- **What this means for the FEP.** The step from "the blanket screens" to "the interior encodes the outside" needs its own condition. This paper is about screening; the encoding question is open beyond these linear checks.
- Limits: linear systems; random networks not designed to encode; linear $R^2$ (exact for Gaussian expectations).

### 3.3 Proposition 2 (exact residual leak for the push score)

**Set-up.** Assume the assumptions of Proposition 1 (block-diagonal positive noise, known start) and a wiring cut, with the push of Eq. 2. Let $\xi = (u, x_I, x_E)$, and set:
- $F = \mathrm{diag}(0, A_{II}, A_{EE})$;
- $Q = \mathrm{diag}(0, D_I, D_E)$;
- $H = [0, A_{BI}, A_{BE}]$.

Given the blanket path, the terms $A_{IB}x_B$ and $A_{EB}x_B$ are known inputs. The filtering covariance then solves
$$
\dot\Sigma = F\Sigma + \Sigma F^{\top} + Q - \Sigma H^{\top} D_B^{-1} H \Sigma, \tag{6}
$$
with initial values below. (The solution exists on $[0, h]$: a filtering Riccati equation with $Q \succeq 0$ and a positive semi-definite start has a global solution.) Initial values:
- $\Sigma_{uu} = I$, $\Sigma_{uI} = \sigma_u I$, $\Sigma_{II} = \sigma_u^2 I$;
- zero for the $E$ rows.

**Reading off the leak.**
- The read is at time $h$, so the filter suffices.
- The numerator of $\ell^c_0(h)$ is $I(u; x_E(h) \mid x(0^-), x_B)$, read off $\Sigma(h)$.
- The cross block between $(u, I)$ and $E$ is driven only by $M$, so $M = 0$ gives a residual leak of exactly zero.
- **When $M \neq 0$, the leak is positive for all $h$ except isolated values.**
  - Near $t = 0$, $\Sigma_{uI} \approx \sigma_u I$ and $\Sigma_{EE} \approx D_E t$, so Eq. 6 gives $\Sigma_{uE}(t) = -\tfrac{t^2}{2}\sigma_u M D_E + O(t^3)$. This is non-zero if and only if $M \neq 0$.
  - $\Sigma(t)$ is real-analytic in $t$, because Eq. 6 has a polynomial right side. So $\Sigma_{uE}$, and with it the numerator of the CPIR, vanishes only at isolated $h$.
  - Conditions: $\sigma_u > 0$; the push acts on all internal coordinates; all of $x_E(h)$ is read; the denominator is positive (Section 2.3). An observer that pushes or reads only some coordinates can see zero leak while $M \neq 0$.
- Its size depends on:
  - $A_{II}$ and $A_{EE}$;
  - $A_{BI}^{\top}D_B^{-1}A_{BI}$ and $A_{BE}^{\top}D_B^{-1}A_{BE}$;
  - $D_I$ and $D_E$;
  - the push, $h$, and the denominator.

### 3.4 Relation to known results
- **Friston et al. (2023), p. 10.** They state the sufficient condition: external states influence only themselves and sensory states, internal states influence only themselves and active states, and external and internal paths are conditionally independent given blanket paths "because there are no flows that depend on both internal and external states (and fluctuations are independent)".
- **Virgo, Rosas & Biehl (2022).** They argue as follows.
  - **The result.** In a time-unrolled sensorimotor loop in discrete time, the current internal and external states are conditionally independent given the history of sensor and actuator states and the initial internal and external states. This follows from d-separation. Their loop has the split structure: sensors are driven by the outside, actuators by the inside.
  - **The start.** They note that the initial states may be correlated, and that for an ergodic system the dependence on them disappears in the infinite-past limit.
  - **Continuous time.** They "expect an analogous result in continuous time".
  - **Relation to this paper.** Proposition 1 gives that continuous-time result for linear diffusions, at the level of whole paths. Beyond it, we add:
    - necessity;
    - the noise weighting: correlated sensory–active noise can break screening even when the influence graph is split (generically it does);
    - non-faithful cancellation;
    - the size of the failure (Proposition 2).
  - **Refinement of their start remark.** With a history before the time asked about, the start effect fades for current states and for paths (Proposition 1b). It stays for paths that begin at the unknown start, however long the path (Section 3.2).
- **Graphical models.**
  - In discrete time, the generic form of Proposition 1 is d-separation on the time-unrolled graph: $b_{t+1}$ is a common child of $s_t$ and $z_t$. $M$ is a noise-weighted, continuous-time form of this collider mechanism (Pearl, 2009; Lauritzen, 1996).
  - "Necessity for weights in general position" is a faithfulness statement (Meek, 1995; Spirtes et al., 2000), and our cancellation case is an unfaithful one.
  - In continuous time, local independence graphs give the graph form of the collider rule (Didelez, 2008; Mogensen et al., 2018; Mogensen & Hansen, 2020) and of correlated noise (Mogensen & Hansen, 2022).
- **Causal graphs for SDEs (Boeken & Mooij, 2026).** They give a d-separation Markov property for sample paths on a finite interval, for additive-noise SDEs with independent initial values (cyclic systems allowed). For a split blanket, Proposition 1's sufficiency is a special case: every inside–outside walk passes the blanket through a non-collider. A blanket state with parents on both sides is a collider that the conditioning opens; $M$ is its noise-weighted, linear form. Their Markov property is one-way (separation implies independence); Proposition 1 adds the converse for linear diffusions, with the exact algebraic condition (including cancellation inside $M$), a start law that need not be independent, and, in Propositions 1b–2, the long-history and causal readers and the size of the leak. We have read their main results, not all of their proofs.
- **Whole time line (Fasen-Hartmann & Schenk, 2025).** For stationary processes, whole component paths are conditionally independent given the rest iff the matching block of the inverse spectral density vanishes for all frequencies. For an Ornstein–Uhlenbeck process, with $S(\omega) = \tfrac{1}{2\pi}(i\omega - A)^{-1}D(i\omega - A)^{-*}$, the $I$–$E$ block of $S^{-1}$ is $2\pi[K + i\omega(C' - G)]$ (the sign of the $i\omega$ term depends on the Fourier convention; the zero condition does not), and for Gaussian processes conditional orthogonality is conditional independence; with $K$, $C'$ and $G$ as in Step 2. With a wiring cut it is $2\pi M$, so for a reader of the whole time line our condition is already known, with its noise weighting. Without a wiring cut, the cancellation $G = C'$, $K = 0$ screens on the whole time line, but not on a finite window, where the endpoint term of Step 4 remains. Example: $A = \begin{pmatrix}-1&0&0.5\\1&-2&1\\0.5&0&-1\end{pmatrix}$, $D = I$, $I = \{1\}$, $B = \{2\}$, $E = \{3\}$: the spectral block is zero at every frequency, while the path dependence on $[0,1]$ with a known start is 0.020 nats (Euler, $dt = 0.02$; post hoc).
- **Filtering.** The Riccati equation, detectability and the forgetting of initial conditions are standard (Kalman & Bucy, 1961; Anderson & Moore, 1979). What we add is their use for path screening (Propositions 1b and 1c).
- **What we add.** For linear diffusions:
  - the finite-window condition (wiring cut and $M = 0$, with the endpoint argument for necessity), complementing the whole-line spectral condition above;
  - the explicit start condition for observers who do not know the start;
  - the long-history and causal-reader results;
  - the exact residual leak.

### 3.5 Numerical checks
**Networks:**
- 8 nodes: $I=\{0,1\}$, sensory $\{2,3\}$, active $\{4,5\}$, $E=\{6,7\}$.
- $I$ and the active states are driven by $I$, sensory and active; sensory and $E$ are driven by $E$, sensory and active.
- Edge probability 0.6, weights U(0.2, 0.8) with random sign; 200 networks.

**Observer:** start known, $\sigma_u = 0.1$, $h = 2$.

| Check | Kind | Prediction | Result |
|---|---|---|---|
| Split blanket: $\ell^c_\delta \propto \delta^p$ | prospectively specified (range from an earlier run on a degenerate case) | $p \in [3.6, 4.4]$ | 3.98 in both directions; observed, not proven |
| Leak vs a both-sides weight $\varepsilon$ | prospectively specified | slope [1.8, 2.2] | 1.99 |
| Leak vs sensory–active noise correlation $\eta$ | prospectively specified | slope [1.8, 2.2] | 1.93 |
| $\varepsilon$ and $\eta$ tuned so that $M = 0$ (direct discrete computation, not the Riccati formula) | algebra check of Prop. 1 | < 5% of the $\eta$-only leak | 0.07%, equal to the no-violation baseline at the same finite $\delta$ |
| $M$ doubled | algebra check | ratio in [3, 5] | 3.98 |
| Riccati, $M = 0$ | algebra check | < $10^{-10}$ | exactly 0 in 200 of 200 |
| Riccati vs direct computation at $\delta = h/32$, $\varepsilon = 0.1$ | algebra check | median deviation < 5% | 0.8% (includes the reference's $\delta^4$ residual) |
| Riccati, leak $\propto \varepsilon^2$ | algebra check | within 5% | 0.19% |

**Notes:**
- The cancellation check uses two different physical violations and an independent computation. It confirms that the zero condition is $M = 0$ and not "no violation". The remaining 0.07% is the finite-$\delta$ residual, also present without any violation.
- The $\delta^4$ rate is a sketch, not verified: the push-driven drift into $B$ between readings has an $O(\delta^2)$ trapezoid error in signal, so $O(\delta^4)$ in information.

### 3.6 How the size depends on system and protocol (prospectively specified)
Same networks; $\varepsilon = 0.1$; Riccati residual leak. Prediction: more than a 2-fold change across $h$ and across $D_E$, and less than a 2-fold change across $\sigma_u$.

| Varied | Group | Values | Median ratio $\ell^c_0$ | Median numerator / denominator (nats) |
|---|---|---|---|---|
| Time scale $h$ | protocol | 0.5, 1, 2, 4 | 0.021, 0.011, 0.0041, 0.0009 (**24×**) | numerator $1.5, 3.2, 2.0, 0.13 \times 10^{-7}$; denominator $1.5, 5.8, 9.3, 3.6 \times 10^{-5}$ |
| External process noise $D_E$ (×) | system | 0.25, 1, 4 | 0.0034, 0.0041, 0.017 (4.9×) | — |
| Push size $\sigma_u$ | protocol | 0.01, 0.1, 1 | 0.0041, 0.0041, 0.0035 (1.2×) | — |

- The ratio column is the median of the per-network ratios. It is not the ratio of the median numerator to the median denominator (for $h = 0.5$ the latter is $0.0099$ against $0.021$; values checked against `M8e.json`).
- Both parts change with $h$. At large $h$, most of the drop comes from the numerator (a 15-fold fall from $h$ = 2 to 4).
- The prediction holds.

## 4. Consequence for DMBD

**Scope.** This section concerns the stated justification of DMBD's path condition, and the reading of its output. It does not concern its inference algorithm, and we did not fit DMBD to data.

### 4.1 What DMBD states and codes
**Statement.** Beck & Ramstead (2025) state $p(s_\tau, z_\tau \mid b_\tau) = p(s_\tau \mid b_\tau)\,p(z_\tau \mid b_\tau)$ (Eq. 7). They present it as a consequence of the dynamics: conditional independence "follows from the absence of any direct causal interaction between s and z".

**Model.** In the paper, boundary dynamics have the form $f_b(s, b, z)$, and the latent noise is diagonal.

**Code.** pyDMBD, commit e058e9c, `models/DynamicMarkovBlanketDiscovery.py`:
- `one_object_mask` (lines 283–302) lets boundary rows depend on $s$, $b$ and $z$, and excludes $s$–$z$ terms;
- the start prior `x0` (lines 66–71) is a full Normal–inverse-Wishart over all latents.

**Inference.** It computes the joint latent posterior exactly (forward–backward), so it does not use Eq. 7.

### 4.2 Check for generic parameters in this model class
**Model:**
- discrete time, $x_{t+1} = Ax_t + w_t$, independent unit noise;
- no $s$–$z$ blocks;
- $s = 2$, $b = 4$, $z = 2$;
- weights U(0.1, 0.4) with random sign.

**The discrete analogue of $M$.** The discrete likelihood has the cross term $s_t^{\top} A_{bs}^{\top}\Sigma_b^{-1}A_{bz} z_t$.

**Test 1** (prospectively specified; $T = 8$, spectral radius ≤ 0.9, 200 models):

| Case | Prediction | Result |
|---|---|---|
| Boundary driven by both sides, start known | > $10^{-3}$ nats in ≥ 95% | 99%; median conditional MI 0.055 nats (22% of the unconditional) |
| Split boundary, start known | < $10^{-10}$ | max $3\times10^{-15}$ |
| Split boundary, stationary start | > $10^{-6}$ in ≥ 90% | 100%; median 0.008 nats (3.3%) |

**Test 2** (prospectively specified; spectral radius ≤ 0.5, 100 models), median MI in nats at $T$ = 8 / 16 / 32:

| Case | Conditional | Unconditional | Verdict |
|---|---|---|---|
| Both sides, start known | 0.038 / 0.082 / 0.170 | 0.135 / 0.319 / 0.687 | Conditional MI grows with $T$ (×4.5; prediction ≥ 3): a steady leak per step, so the failure is structural |
| Split, stationary start | 0.0043 / 0.0043 / 0.0043 | 0.163 / 0.329 / 0.652 | Constant in $T$ (×1.00; prediction ≤ 1.5): a start effect only |
| Split, start known | < $4\times10^{-14}$ | 0.120 / 0.274 / 0.583 | Exact screening |

**Two samples.** The split, stationary-start median here (0.0043) and in Sections 3.2 and 3.2b (0.0048) come from two different samples of 100 models from the same class and seed stream: M10 drew 100 models of the first case before the split models.

**On the thresholds.** With random weights, exact factorisation fails except on a measure-zero set where the discrete $M$ cancels. The thresholds therefore report size, not existence.

### 4.3 Reading
- **Necessary, not sufficient.** No direct interaction is necessary for Eq. 7, but not sufficient. Also needed:
  1. a discrete $M = 0$, for example a split boundary;
  2. a start law that factorises given $b_0$.
- **The failure grows.** For generic parameters in DMBD's model class, the failure is structural and grows with the length of the path.
- **What is affected.** DMBD's inference is unaffected. The reading of a detected partition as a path-level blanket is affected.
- **Diagnostic.** Report the discrete-time $\hat M = \hat A_{bs}^{\top}\hat\Sigma_b^{-1}\hat A_{bz}$, and the leak computed from the discrete-time conditional covariance (not the continuous-time Riccati formula).
- **Remedy.** Enforce a split mask on the boundary rows and a start prior that factorises given $b_0$.

## 5. What the observation protocol adds

### 5.1 Observers who do not read the blanket: leak grows with time scale

**Proposition 3.** Assume a wiring cut with $\|A_{BI}\|_F > 0$. Then
$$
\ell_{I\to E}(h) = \frac{h^2}{4}\kappa_{IE} + O(h^3), \qquad \kappa_{IE} = \frac{\|A_{EB}A_{BI}\|_F^2}{\|A_{BI}\|_F^2} \le \sigma_{\max}(A_{EB})^2 . \tag{7}
$$
More generally, let $d_E$ and $d_B$ be the shortest path lengths from $I$ to $E$ and to $B$, and assume the leading matrix-power terms are not cancelled. Then the leak grows as $h^{2(d_E-d_B)}$. For example, a blanket of two layers, $I \to B_1 \to B_2 \to E$, gives $h^4$.

This is a statement about the total-effect score, not about the blanket itself.

| Algebra check (200 networks) | Result |
|---|---|
| Slope of the leak in $\log h$, one layer | 2.000 |
| Slope of the leak in $\log h$, two layers | 4.000 |
| Error of the leading term at $h = 0.1$ | 0.9% |
| Blanket time (leak 0.5) predicted as $2/\sqrt\kappa$, from equating the leading terms of numerator and denominator of the share (the leading term alone, $\tfrac{h^2}{4}\kappa = 0.5$, would give $\sqrt{2/\kappa}$) | median prediction/exact 0.97 |

These are blanket analogues of time-scale effects in undersampled causal discovery (Danks & Plis, 2013; Barnett & Seth, 2017).

### 5.2 Proposition 4 (slow limit)

**Assumptions:**
- the eigenvalue of $A$ with the largest real part is real, simple and strictly dominant, with right eigenvector $v$ and left eigenvector $w$;
- $w_I \neq 0$ and $w_E \neq 0$;
- $v_B \neq 0$.

**Claim.** (The formula holds in whatever units $A$ is written; the score $\rho$ is what depends on units.) In noise units,
$$
\rho(h\to\infty) = \frac{\|v_B\|^2}{\|v_B\|^2 + \max(\|v_I\|^2,\|v_E\|^2)} . \tag{8}
$$

**Complex dominant pair** $\alpha \pm i\omega$ (strictly dominant):
- remove the common envelope, $\tilde R(h) = e^{-\alpha h}R(h)$;
- average the numerator and the denominator of the share separately over one period, as $h \to \infty$;
- take the ratio.

This defines a variant of the score, not its pointwise limit. Eq. 8 then holds with $|v|$. *Proof.* The pair contributes $2e^{\alpha h}\mathrm{Re}(v w^{*} e^{i\omega h})$ to $R(h)$, so each entry is $2e^{\alpha h}|v_a||w_j|\cos(\omega h + \phi_a - \psi_j)$. Averaged over one period, its square is $2e^{2\alpha h}|v_a|^2|w_j|^2$: the phases drop out entry by entry. The averaged sums over blocks therefore factorise as $\|v_E\|^2\|w_I\|^2$ and $\|v_B\|^2\|w_I\|^2$, and $\|w_I\|^2$ cancels as in the real case. $\square$ (Algebra check: max deviation $4\times10^{-8}$ in 23 networks.)

**Reading.**
- A slow total-effect observer ranks partitions by where the slowest mode lives, not by the wiring.
- The sparse-coupling and slow-mode readings of blankets are the fast and slow limits of one family of observers.
- This is the blanket analogue of Markov stability for communities (Delvenne et al., 2010; Lambiotte et al., 2014), and it connects to the adiabatic reduction of Friston (2019).

| Algebra check (12-node networks, 24 with a real dominant mode) | Result |
|---|---|
| Rank correlation of $\rho(h)$ with $\rho(\infty)$ at $h$ = 1, 3, 10, 30 | 0.38, 0.77, 0.998, 1.000 |

### 5.3 Does the best partition depend on the time scale?
The computations are exact. We report how often the best (argmax) partition changes between time scales, and compare large with small changes of $h$, as a dose–response check on how sensitive the argmax is.

**Total-effect observers** (prospectively specified; 50 random 12-node networks; all 4/4/4 partitions).

| Representation | Best differs, $h$ = 0.5 vs 3 | Best differs, $h$ vs $1.2h$ | Ratio |
|---|---|---|---|
| Noise units (50 networks) | 98% | 40% | 2.4 |
| Random units (3 per network, 150) | 98% | 45% | 2.2 |

- In an earlier run, the direct couplings ($h \to 0$) and $h = 3$ gave different best partitions in 99.5% of 200 networks.
- The prediction that the large change gives at least 3 times the small one was not reached in either representation.
- **Reading.** The total-effect argmax is very sensitive to the time scale; even a 20% change moves it in 40–45% of networks. The pattern is the same in random units, so it is not a unit artefact.

**Observers who read the blanket** (prospectively specified).
- **Design:**
  - read interval $\delta = 0.1$; time scale 0.5 vs 3; the CPIR for both observers;
  - "differs" means that the argmax among eligible candidates changes; eligible candidates are those in the top 10% under both observers;
  - at least 30 eligible networks required.
- **Networks:** 150 split-blanket networks (8 nodes, 210 candidates each), each with weak violations of one type:
  - "direct": internal–external edges;
  - "split": a blanket state hears the other side.

| Networks | Differs (0.5 vs 3) | Differs ($h$ vs $1.2h$) | Paired difference, 95% CI | Planted partition best (0.5 / 3) |
|---|---|---|---|---|
| All (150) | 29.3% | 11.3% | +18 points (11–25) | 52% / 51% |
| Direct (75) | 42.7% | 16.0% | +27 (15–37) | 23% / 32% |
| Split (75) | 16.0% | 6.7% | +9 (1–19) | 81% / 71% |

**Verdict.**
- **Dose–response is present:** a large change of time scale moves the argmax more often than a small one (the paired CI excludes zero).
- The prospectively specified size (at least 30% and at least 3 times) was not reached.
- "Direct > split" holds. The larger effect with direct edges is partly a normalisation effect (Section 2.3).
- With $M$-type imperfections only, observers mostly agree on the planted partition.

### 5.4 Critical slowing: where the slow mode lives matters

**Set-up:**
- split-blanket networks with no violations, so $M = 0$;
- a positive two-cycle is raised until the slowest eigenvalue reaches $-0.02$;
- **main case:** the loop is inside $I$;
- **location control** (prospectively specified, same slowness): the loop is inside the sensory pair.

| Loop | Networks | Share of slow mode on $I$ / $B$ | Median change of total-effect $\rho(3)$ per network |
|---|---|---|---|
| Inside $I$ | 87 | 0.72 / 0.20 | −0.19 |
| Inside $B$ (sensory pair) | 92 | 0.16 / 0.62 | +0.07 |

**Predictions:**

| Prediction | Result | Verdict |
|---|---|---|
| Inside-$I$ loop: drop ≥ 0.2 | 0.19 per network (0.25 as a difference of medians in the first run; Appendix A) | not reached |
| Inside-$I$ loop: share on $I$ ≥ 0.8 | 0.72 | not reached |
| Inside-$B$ loop: no drop (≥ −0.02) | +0.07 | holds |
| Inside-$B$ loop: share on $B$ ≥ 0.6 | 0.62 | holds |
| Blanket reader every $\delta = 0.1$: $\rho^c$ | $1 - 2\times10^{-6}$; the finite-read residual (falls as $\delta^4$); the continuous limit is exactly 1 by Prop. 1 | algebra check |

**Reading.**
- Slowing the inside makes the same blanket fade for a slow total-effect observer.
- Slowing the blanket makes it look better.
- The location of the slow mode decides the reading, as Proposition 4 says.

### 5.5 Conditional-independence blankets and push observers near equilibrium

**Definitions.** A **snapshot blanket** is a zero block in the stationary precision: $I$ and $E$ are conditionally independent given $B$ in the stationary density.

**Set-up.**
- Take $A = -(D + sQ)\Pi$, with $D$ diagonal and $Q$ antisymmetric, and let the **diffusion matrix be the same $D$** (Eq. 1).
- Then $A\Sigma + \Sigma A^{\top} = -D$ holds for $\Sigma = \tfrac12\Pi^{-1}$. So the stationary precision is $2\Pi$, and its zero pattern is that of $\Pi$.
- If the diffusion matrix differs from $D$, this fails. Check M2b (200 cases): with diffusion $D$, the precision equals $2\Pi$ to $6\times10^{-15}$. With diffusion $I$ while $D \neq I$, it differs in 100% of cases, and its $I$–$E$ block is non-zero.
- Let $\Pi_{IE} = 0$.

**Result.** A fast push observer sees a direct leak
$$
s^2\|(Q\Pi)_{EI}\|^2/\|(D\Pi)_{BI}\|^2 + O(s^3).
$$
- It is zero at equilibrium ($s = 0$) and generically positive otherwise; $(Q\Pi)_{EI} = 0$ is possible for $s \neq 0$. It needs $(D\Pi)_{BI} \neq 0$, and $D$ diagonal so that $(D\Pi)_{EI} = 0$.
- Algebra check: slope 1.99; formula/exact = 1.001 at $s = 0.005$.
- Aguilera et al. (2022, Eq. 29) start from a sparse drift and expand the precision. We start from a sparse precision and expand the drift.

### 5.6 Layered systems: thickness gives the split; several nested exact shells; a moving blanket

**The system.**
- A layered ("onion") system: a core $L_0$, then layers $L_1, \dots, L_4$.
- Each state interacts only within its layer and with the adjacent layers.
- Diagonal noise, so $D_B$ is diagonal.

**Lemma (thickness gives the split).** Take the inside to be the inner layers and the outside to be the outer layers.
1. **One layer.** Suppose some state of the blanket $B = L_k$ has edges to both sides. Then $M \neq 0$ for generic weights.
2. **Two layers.** For $B = L_k \cup L_{k+1}$ there is a wiring cut, and $M = 0$ for **all** weights. The inner layer hears only the inside and $B$; the outer layer hears only the outside and $B$.

**Scope.** The Lemma is about nearest-neighbour layered graphs with diagonal noise. In general graphs, a one-layer blanket can have $M = 0$ (sparse or cancelling weights, or separate sensory-like and active-like states within one layer), and a two-layer shell can have $M \neq 0$ (links that skip a layer, or correlated noise). Calling the two faces "sensory" and "active" is a reading of the Lemma, not part of it.

**Who sees an exact blanket.** By Propositions 1, 1b and the start law (Section 3.2), a two-layer shell is exact for continuous readers who either:
- know the start; or
- have a start law that factorises given $x_B(0)$; or
- are in the long-history limit (Proposition 1b; exponentially close at a finite history).

A continuous reader who starts at a stationary state without that history does see a small start leak.

**Position selects among exact shells.** Consider observers that push and read one layer and read the candidate blanket continuously, in the class above.
- **Exactness does not depend on where the observer stands.** Every two-layer (or thicker) shell is exact.
- **Position decides which shells lie between the observer and the core:**
  - a near observer (access to $L_3$) meets one two-layer shell, $L_1 \cup L_2$;
  - a far observer (access to $L_4$) meets two nested ones, $L_1 \cup L_2$ and $L_2 \cup L_3$.
- **What this shows and does not show.**
  - It shows that several valid partitions exist at once.
  - It gives no principle that selects one of them as "the" particle (open).

**A moving blanket.** Let the core grow a direct path $L_0 \to L_2$ of strength $\varepsilon$.
- The inner shell $L_1 \cup L_2$ now has a layer that hears both sides, so $M \neq 0$ and it leaks.
- The outer shell $L_2 \cup L_3$ stays exact.

**Checks** (M15; 50 networks; start known; prospectively specified):

| Check | Prediction | Result |
|---|---|---|
| One-layer blanket leaks | ≥ 95% | 96% (median leak 0.13) |
| Two-layer shells, $M$ and leak | 0 in all | 0 and 0 |
| Near and far observers' exact shells | as stated | as stated |
| Inner-shell leak at $\varepsilon$ = 0, 0.1, 0.2, 0.3 | rises; slope in $\log\varepsilon$ in [1.6, 2.4] | 0, 0.041, 0.14, 0.29; slope 1.77 |
| Outer-shell leak under the shortcut | 0 | 0 |

One network without a causal route from the core to the outside was excluded. Its score is undefined; the code was fixed after the first run to report such cases.

### 5.7 A real case: fMRI effective connectivity

**Data and fit.**
- Resting and movie-viewing fMRI from Gilson et al. (2018): 22 subjects, 5 runs each, 66 regions, 300 volumes per run, TR 2 s.
- We fitted multivariate Ornstein–Uhlenbeck models with the authors' code (pyMOU) and settings.

**Partitions** (fixed in advance):
- **Inside:** left or right visual cortex, or left or right sensorimotor cortex.
- **Blanket:** all regions connected to the inside.
- **Outside:** the rest.
- The fitted drift has no direct inside–outside term, so the wiring criterion accepts each partition (a statement about sparsity, not about screening).

**Score and null.**
- The score is the continuous reader's residual leak ($h$ = 2 TR).
- The null comes from the same pipeline applied to data simulated from each fitted model, with the blanket split ($M = 0$ by construction).
  - **Split model:** each blanket state keeps inputs only from the side (inside or outside) with the larger summed $|A|$; the fitted noise is kept. If a split model was not stable, its diagonal was shifted so that its largest eigenvalue real part was $-0.01$ (the same rule in every case; it changes the null's time scales).
  - **One null per real case:** one 300-volume series is simulated (Euler, $dt$ = 0.02 TR, burn-in 200 TR) and passed through the same filter, normalisation and fit.
  - **Cases:** 440 = 22 subjects × 5 runs × 4 inside sets. They are not independent (subjects, runs and overlapping partitions); we did no subject-level analysis.

| | Real data | Null ($M = 0$) |
|---|---|---|
| Median leak | 0.33 | 0.51 |
| Real leak > null leak | 18% of 440 cases (prediction ≥ 80%) | — |

**The prospectively specified failure condition fired.** With 300 volumes, this pipeline cannot tell the estimated leak apart from its null.

**Post-hoc simulation** (10 fitted models as ground truth):

| Volumes per fit | True-model leak ranked above split-model leak | Median refitted leak, true / split (true values 0.26 / 0) |
|---|---|---|
| 300 | 52% | 0.55 / 0.54 |
| 1200 | 88% | 0.58 / 0.35 |
| 4800 | 90% | 0.56 / 0.28 |

**Reading.**
- **No decision.** The present estimator and calibration do not support a decision about screening in these recordings. Estimation bias is large: the refitted leak of a model with true leak 0 stays well above 0 (0.54 at 300 volumes). We have not excluded other causes of the real-vs-null difference: null dynamics that do not match the data, preprocessing, regularisation, or model misspecification.
- **What is tested.** Only screening in the fitted Ornstein–Uhlenbeck model. Screening in the measured BOLD process, or in the underlying neural dynamics, is a further step; band-pass filtering and 2-s sampling are far from continuous reading.
- **Not a data requirement.** Longer series separate the two model classes better in this simulation. This is not a general data requirement for brain recordings.
- **Two different tasks.**
  - **Rejecting an exact zero:** showing $M \neq 0$.
  - **Supporting an approximate zero:** showing that the leak is below a stated margin.
  - The second needs valid uncertainty and error control relative to a margin fixed in advance (an equivalence test). A bias-corrected estimator is one route. Neither is available yet.

## 6. Discussion

### 6.1 What can be agreed on
Fix a representation (noise units, or information-based scores) and the protocol class of Proposition 1 (start known, continuous reading of the blanket).
- **Exactness.** All observers in this class agree on whether screening is exact: a wiring cut and $M = 0$.
- **Finite read intervals.** In our examples, a reader at a finite interval saw a small leak even when $M = 0$, and it fell as $\delta^4$ (observed, not proven). Systems with special structure can stay screened at finite intervals.
- **Limit.** What such observers can agree on is whether this leak vanishes in the continuous-read limit.
- **Size and ranking.** Observers need not agree on the size of a failure (Section 3.6). Observers who read the blanket at a fixed finite interval ($\delta = 0.1$, so outside the continuous-read class) agree on the ranking of partitions to a large extent when the imperfections are of the $M$ type (Section 5.3).

We suggest reporting a blanket together with the protocol class for which it is claimed. This gives a measurable version of the question whether blankets are modelling tools or real boundaries (Bruineberg et al., 2022).

### 6.2 Limitations
- Linear Gaussian dynamics only; in nonlinear systems $M$ may depend on the state.
- Small networks with random weights.
- The noise is block-diagonal across groups.
- The $\delta^4$ rate is observed, not proven.
- Several prospectively specified sizes were not reached (Sections 5.3 and 5.4).
- DMBD was not fitted to data.
- One real-data case only, with no decision; no calibrated test of the leak yet.
- Screening only. Whether the internal state encodes the outside is a separate question; the synchronisation map needs its own condition; for the causal reader no exact map was found in 90 networks, while for the snapshot reading it held in the small-blanket controls (Section 3.2d).
- No selection principle among several exact shells (Section 5.6).

### 6.3 Conclusion
- **Start known.** For linear diffusions with block-independent, positive-definite noise and a known start, finite-window path screening by a continuously read blanket holds if and only if there is a wiring cut and $M = A_{BI}^{\top}D_B^{-1}A_{BE} = 0$.
- **Start unknown.** Under the same assumptions, a start law (with a positive density) that factorises given $x_B(0)$ is also needed; given a wiring cut and $M = 0$ it is necessary and sufficient. With a long history of blanket readings, it is not needed (Proposition 1b).
- **Causal reader.** A reader of the past blanket only needs the same $M = 0$ for the current states (Proposition 1c).
- **Encoding is separate.** A synchronisation map needs its own condition (Section 3.2d).
- **A split blanket** guarantees $M = 0$.
- **When screening fails**, the residual leak can be computed exactly, and its measured size depends on system and protocol together.
- **For detection methods.** Methods that state a path-level blanket condition should check $M$ and the start law, and report the protocol class for which a partition is claimed.

## Appendix A. Disclosure: changes, failures and post-hoc analyses
1. **Collider diagnosis.** An earlier claim, that continuous reading screens off every wiring cut, failed its prospectively specified test (the leak stayed at 0.12). The collider diagnosis was post hoc; it was then checked prospectively (Section 3.5).
2. **Necessity argument.**
   - v0.2 used the push observer's start covariance. v0.4 gives the path-likelihood argument and the no-push filtering term.
   - The filtering term's check missed its 1% tolerance at $t = 0.01$ (2.0%). A post-hoc check showed the expected linear shrinkage.
3. **Corrected statements from earlier drafts:**
   - "the leak depends on the system only through $M$";
   - "with diagonal noise $M = 0$ means a split";
   - "a split gives $M = 0$" (without requiring uncorrelated sensory–active noise);
   - "the system decides whether" (stated without the protocol class);
   - $D_E$ filed under "observer";
   - the floor named as if it were a lower bound;
   - a units check that could not discriminate (v0.3; redone in Section 5.3);
   - the time-scale comparison called a "null" (now a dose–response check).
4. **Time-scale tests.**
   - An early test confounded time scale with read interval.
   - Another had too few eligible networks (4 of 40), and its minimax prediction failed.
   - The final test changed the small step from $1.1h$ to $1.2h$ before its run.
   - The paired bootstrap is post hoc.
5. **Critical slowing.**
   - The first null loop (a negated symmetric loop) was itself a positive cycle; it was replaced before the main measure was read.
   - The first run summarised the drop as a difference of medians (0.25); the location run used the median change per network (0.19). The verdict "not reached" uses the second.
6. **Earlier analyses.** A check on blanket thickness (Section 5.1) excluded networks without a length-3 path after a first run. Its blanket-time threshold of 0.5 was untestable in about half of the networks.
7. **Added after internal reviews:**
   - DMBD scale and length analyses;
   - leak-size sensitivity and nats;
   - the units dose–response check;
   - the location control.
   These were prospectively specified where stated.
8. **Virgo et al. (2022)** was read in full for v0.5. Section 3.4 now states their result exactly. The current-state vs path check (Section 3.2) was prospectively specified for $T$ = 8, 16, 32. At $T = 8$ the current-state dependence was already about 0, so its ratio could not be evaluated, and the short histories ($T$ = 1, 2, 4) were added post hoc.
9. **Changes after the Emmy team review of v0.5 (v0.6).**
   - Necessity is now proved in one step for both the wiring cut and $M$, with the endpoint term (the earlier Step 4 did not exclude cancellation).
   - Novikov's condition was replaced by Beneš's condition.
   - The noise notation of Section 5.5 is corrected.
   - Proposition 1b is new, with a full proof.
   - Sections 5.6 and 5.7 now carry the layered-system and real-data material that was previously only in the companion paper.
10. **Failed prediction (M13b c).** We predicted that one interior mode hidden from the blanket would slow forgetting. It did not. The explanation (the sum of the two sides' rates) was derived after the run. It was then checked post hoc (M13c): hidden slow modes on both sides do slow forgetting.
11. **Real data (Section 5.7).** The prospectively specified test failed (18% vs ≥ 80%). The data-length simulation is post hoc and covers only 10 models.
12. **Layered systems (Section 5.6).** After the first run, the code was fixed to report networks with no causal route, instead of letting them corrupt the medians.
13. **Changes after the second review round (v0.7; Vector, Nova and Trace on the main paper v0.6, with this note).**
   - Proposition 1b, step (iv): the data-processing argument was wrong; it is replaced by the chain rule, which gives equality with the start dependence. The rate is stated as an asymptotic bound; "whatever the stationary density" is restricted to the Gaussian stationary start.
   - Proposition 1: the necessity of the wiring cut and $M = 0$ is extended to an unknown start; two rigour remarks added (Wiener integrals; the almost-sure shift identity).
   - Proposition 2: a proof that the leak is positive for all but isolated $h$ when $M \neq 0$, with the conditions on the push and the read.
   - Proposition 1c (causal reader) is new, with check M17; Section 3.2d (encoding) is new, with check M17b.
   - The two samples behind 0.0043 and 0.0048 are explained; current states vs paths is stated with the intervals; Section 5.5 says "generically"; Section 5.6 states its scope; Section 5.7 documents the null and no longer says that estimation bias is the cause.
14. **M17 and M17b.** M17 was prospectively specified; its three predictions passed. Its encoding measure (internal state) was exploratory. Trace's review judged it the important result; we then wrote M17b, a fairer test of the FEP's actual claim (internal expectation, synchronisation map), with a positive control, before running it. All three M17b predictions passed.
15. **Changes after the third review round (v0.8; Vector, Nova on the main paper v0.7; Trace on this note v0.7).** Trace found no blockers. Stated assumptions added (positive start density; $\Sigma_{BB} \succ 0$; infinite-past limit; wiring cut assumed for the causal reader); the $M \neq 0$ case of Proposition 1b is now proved; "hidden" modes are now "slow closed-loop modes"; unit invariance is stated for a push transformed with the variables; Section 5.4 labels $\rho^c$ correctly.
16. **M17c.** Prospectively specified. Part (a) passed. Part (b) had two code problems: in run 1 the generator did not force the collider (some networks had $M = 0$), and run 2 stopped on a rank-deficient $W_I$. After both fixes (disclosed in the script; run-1 output kept), (b) gave 29 of 30, so its prediction "in all 30" failed; the one exception has rank-deficient $W_I$, which the kernel condition predicts. The tool now uses a pseudo-inverse in the encoding score.
17. **Changes after the fourth review round (v0.9).** Trace: ready, no blockers; Vector: minor; Nova: major on framing. Changes: density condition in Section 3.2; the active-state finding scoped and renamed (proprioception-like); the encoding non-implication labelled by reading; the Section 5.4 residual explained as a finite-read residual (it was mislabelled as the continuous reader's); appendix numbering fixed.
18. **M17d.** Post hoc and exploratory, after the fourth review: causal-reader encoding in 90 networks (no exact map found).
19. **Release candidate (v0.15).** A prior-art search found Fasen-Hartmann & Schenk (2025), which gives the whole-time-line spectral form of the screening condition; Section 3.4 now cites it and narrows our claim to the finite window. The example in Section 3.4 (M16b) is post hoc. Da Costa et al. (2021, Lemma 2.1) is credited for the kernel condition. References were checked; the remaining * items are marked.
20. **Removed material.** Push-and-sham calibration and wind-tunnel tests from earlier drafts are not part of this paper. They will be reported separately.

## Reproducibility
All paths are relative to the code archive `fep-blanket-code-v1.0.zip`, published together with this note (also at github.com/casualwriter/ifo-fep). Each script prints its results and writes one JSON file. Predictions are in the script headers.

| Script (`code/track-a/`) | Section |
|---|---|
| `M3_loop2.py`, `M3b_collider_posthoc.py` | 2.3–2.4, App. A |
| `M4_loop3.py`, `M4b_cancel.py` | 3.5 |
| `M8_riccati_floor.py`, `M8b_sensitivity.py`, `M8c_nats_necessity.py`, `M8d_necessity_order.py`, `M8e_nats_by_h.py` | 3.1, 3.3, 3.5, 3.6 |
| `M9_dmbd_path_factorisation.py`, `M9b_scale.py`, `M10_dmbd_transient.py`, `M10b_uncond.py`, `M10c_current_vs_path.py`, `M10d_short_T.py` | 3.2, 4 |
| `M1_math_checks.py`, `M2_loop1.py` | 5.1, 5.2, 5.5 |
| `T4_observer_dependence.py`, `M11_units_robust.py`, `M11b_units_null.py` | 5.3 |
| `M7_fair_h_test.py` (paired bootstrap in `M7b_paired.json`) | 5.3 |
| `M6_loop5.py`, `M12_loop_location.py` | 5.4 |
| `M16_endpoint.py`, `M16b_whole_line.py` | 3.1, 3.4 |
| `M13_long_history.py`, `M13b_detectability.py`, `M13c_two_sided.py`, `M13d_late_slope.py` | 3.2b |
| `M17_causal_reader.py`, `M17b_sync_map.py`, `M17c_sensory_only.py`, `M17d_causal_encoding.py` | 3.2c, 3.2d |
| `code/tools/blanket_report.py` (functions `causal_leak`, `sync_r2`) | 3.2c, 3.2d |
| `M2b_noise_notation.py` | 5.5 |
| `M15_layers.py` | 5.6 |
| `M14_real_mouec.py`, `M14b_power.py` (data: the public copy at github.com/mb-BCA/notebooks_review2019, Apache-2.0; pyMOU commit f0626be; neither is included in the archive) | 5.7 |

Python 3.13, NumPy, SciPy. Each run takes under 25 minutes on an old laptop. DMBD code: `code/dmbd/pyDMBD/`, commit e058e9c, unchanged.

## References
Items marked * have details we could not fully confirm (a corollary number; a workshop paper).

- Anderson, B.D.O. & Moore, J.B. (1979). *Optimal Filtering*. Prentice-Hall.
- Aguilera, M., Millidge, B., Tschantz, A. & Buckley, C.L. (2022). How particular is the physics of the free energy principle? *Physics of Life Reviews* 40, 24–50. doi:10.1016/j.plrev.2021.11.001
- Barnett, L. & Seth, A.K. (2017). Detectability of Granger causality for subsampled continuous-time neurophysiological processes. *Journal of Neuroscience Methods* 275, 93–121.
- Beck, J. & Ramstead, M.J.D. (2025). Dynamic Markov blanket detection for macroscopic physics discovery. arXiv:2502.21217. Code: github.com/bayesianempirimancer/pyDMBD.
- Boeken, P. & Mooij, J.M. (2026). Causal graphs, Markov properties and do-calculus for stochastic differential equations. arXiv:2607.12140.
- Bruineberg, J., Dołęga, K., Dewhurst, J. & Baltieri, M. (2022). The Emperor's New Markov Blankets. *Behavioral and Brain Sciences* 45, e183.
- Da Costa, L., Friston, K., Heins, C. & Pavliotis, G.A. (2021). Bayesian mechanics for stationary processes. *Proceedings of the Royal Society A* 477(2256), 20210518. doi:10.1098/rspa.2021.0518
- Danks, D. & Plis, S. (2013). Learning causal structure from undersampled time series. NIPS 2013 Workshop on Causality (workshop paper).*
- Delvenne, J.-C., Yaliraki, S.N. & Barahona, M. (2010). Stability of graph communities across time scales. *PNAS* 107(29), 12755–12760.
- Didelez, V. (2008). Graphical models for marked point processes based on local independence. *Journal of the Royal Statistical Society B* 70(1), 245–264. doi:10.1111/j.1467-9868.2007.00634.x
- Friston, K. (2019). A free energy principle for a particular physics. arXiv:1906.10184.
- Fasen-Hartmann, V. & Schenk, L. (2025). Partial correlation graphs for continuous-parameter time series. *Metrika* 88(6), 1425–1460. doi:10.1007/s00184-025-00992-6
- Friston, K., Da Costa, L., Sakthivadivel, D.A.R., Heins, C., Pavliotis, G.A., Ramstead, M. & Parr, T. (2023). Path integrals, particular kinds, and strange things. *Physics of Life Reviews* 47, 35–62. arXiv:2210.12761v3.
- Gilson, M., Deco, G., Friston, K.J., Hagmann, P., Mantini, D., Betti, V., Romani, G.L. & Corbetta, M. (2018). Effective connectivity inferred from fMRI transition dynamics during movie viewing points to a balanced reconfiguration of cortical interactions. *NeuroImage* 180, 534–546.
- Heins, C. & Da Costa, L. (2022). Sparse coupling and Markov blankets: a comment on "How particular is the physics of the free energy principle?" by Aguilera, Millidge, Tschantz and Buckley. *Physics of Life Reviews* 42, 33–39. doi:10.1016/j.plrev.2022.06.001
- Kalman, R.E. & Bucy, R.S. (1961). New results in linear filtering and prediction theory. *Journal of Basic Engineering* 83(1), 95–108.
- Karatzas, I. & Shreve, S.E. (1991). *Brownian Motion and Stochastic Calculus* (2nd ed.). Springer.*
- Lambiotte, R., Delvenne, J.-C. & Barahona, M. (2014). Random walks, Markov processes and the multiscale modular organization of complex networks. *IEEE Transactions on Network Science and Engineering* 1(2), 76–90.
- Lauritzen, S.L. (1996). *Graphical Models*. Oxford University Press.
- Meek, C. (1995). Strong completeness and faithfulness in Bayesian networks. *Proceedings of UAI 1995*, 411–418.
- Mogensen, S.W., Malinsky, D. & Hansen, N.R. (2018). Causal learning for partially observed stochastic dynamical systems. *Proceedings of UAI 2018*, PMLR, 350–360.
- Mogensen, S.W. & Hansen, N.R. (2020). Markov equivalence of marginalized local independence graphs. *Annals of Statistics* 48(1), 539–559. doi:10.1214/19-AOS1821
- Mogensen, S.W. & Hansen, N.R. (2022). Graphical modeling of stochastic processes driven by correlated noise. *Bernoulli* 28(4), 3023–3050.
- Pearl, J. (2009). *Causality* (2nd ed.). Cambridge University Press.
- Spirtes, P., Glymour, C. & Scheines, R. (2000). *Causation, Prediction, and Search* (2nd ed.). MIT Press.
- Virgo, N., Rosas, F.E. & Biehl, M. (2022). Embracing sensorimotor history: Time-synchronous and time-unrolled Markov blankets in the free-energy principle. *Behavioral and Brain Sciences* 45, e215. doi:10.1017/S0140525X22000334 (accepted manuscript read: University of Sussex repository)

## Version history
- **v1.0.1 (October 2026).** Adds Boeken & Mooij (2026) and the comparison with their d-separation Markov property for SDEs (an external AI review, Stanford Agentic Reviewer, pointed to it) (Section 3.4). The main paper 1.0.1 also states the A3 proof sketch in one place and adds a units recommendation. No result changes.
- **v1.0 (October 2026).** First public release (Zenodo). Content as working version v0.15, with the corrections in the last item.
- Working versions v0.1–v0.15 (July–October 2026) were revised after reviews by the Emmy team and by AI reviewers. The changes, failed predictions and corrected claims are listed in Appendix A.
- Before release, version 1.0 was corrected after a last review by two AI reviewers (Nova and Vector, 20 October 2026; both: minor revision): the companion's title no longer says that a split is needed; its conclusion keeps the wiring cut in the iff; claims about many blankets are limited to structured systems; diagnostics are kept apart from screening; the tool's defaults are stated. A final check by Trace found nothing wrong (ready for release); its wording points were applied.
