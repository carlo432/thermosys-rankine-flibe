# Verification, validation, and what's actually been established

Three different claims get conflated in "is the physics right?". They are separated here,
because only two of them are settled.

| Claim | Question | Status |
|---|---|---|
| **Verification** | Are the equations solved correctly? | **Established** — three independent ways |
| **Validation** | Are they the right equations? | **Partial** — strong for water, provenance-only for FLiBe |
| **Uncertainty** | How much does the answer move if the inputs are wrong? | **Quantified** — see §3 |

---

## 1. Verification

### 1.1 Internal — falsifiability (`FALSIFIABILITY.md`)

34 physics terms, each deleted in turn by `tools/mutation_check.py`, each killed by a named
test. This proves the tests can fail. It does
**not** prove any number is right.

Two mutations initially survived, and both were real gaps, not curiosities:
- the loop's convergence residual guard — never exercised because the solver always converged
- the boiler's zone-clipping branch — only reachable at low heat input, which no test used

Both now have targeted tests. That is the harness doing its job.

### 1.2 Independent formulation — IAPWS-IF97 (`validation/cross_check_if97.py`)

Every cycle recomputed on a different property formulation (IF97 industrial vs CoolProp's
IF95 Helmholtz) with the cycle algebra written a second time from scratch, importing nothing
from `thermosys`.

| Anchor | Worst deviation |
|---|---|
| 5 saturation points (T_sat, h_f, h_g) | 1.9e-4 |
| 10 cycle results (η, x, y across all four cycles) | 3.1e-4 |

IF97 and IF95 agree to ~1e-4–1e-3 by construction, so these deviations are the two
formulations differing, not an error in either. A flipped isentropic relation or a wrong
state point would show up orders of magnitude larger. **This verifies the cycle
construction; it cannot confirm the anchors at the 1e-6 tolerance the suite asserts.**

### 1.3 Third-party published values (`tests/test_external_references.py`)

The only checks here that compare against numbers *someone else* published:

| Case | Published | Computed |
|---|---|---|
| Reheat 8 MPa → 4 MPa → 10 kPa, 500 °C: HP exit | 3185 kJ/kg | 3184.6 |
| same, LP exit | 2247 kJ/kg | 2246.8 |
| Ideal Rankine 15 MPa/600 °C/10 kPa: exit quality | 80 % | 0.8041 |
| same with reheat at 1 MPa: exit quality | 98 % | 0.9843 |

Worth being precise about a limitation: my own headline cases (3 MPa saturated → 75 kPa;
8 MPa/500 °C with 1 MPa reheat; regeneration at 1 MPa extraction) are *not* digit-for-digit
published examples — the standard textbook versions use different conditions. So these four
verify the cycle machinery at nearby conditions rather than certifying my exact anchors.
A remaining independent check is comparison against EES state points for the same conditions.

---

## 2. Validation — where the two fluids genuinely differ

**Water: validated by construction.** IAPWS-95/IF97 are fits to experimental data, and the
cycle results agree with two independent implementations of them. There is no meaningful
open question on the steam side.

**FLiBe: provenance only, and extrapolated.** The correlations come from Sohal/INL and
Attarian, which are themselves fits to measurements — but:

- viscosity is fitted over **873–1073 K** and used over **732–1703 K**
- conductivity is fitted over **800–1200 K** and used to 1703 K
- no loop-level experimental comparison exists in this repo

No amount of additional testing fixes that; it is a data problem, not a code problem. Real
validation would need molten-salt loop data (MSRE-class experiments). What *can* be done
without new data is quantify the consequence, which is §3.

---

## 3. Uncertainty propagation (`thermosys/analysis/uncertainty.py`)

The loop's `UA` is computed from the coolant's own properties via Dittus-Boelter, so
viscosity and conductivity actually enter the answer. Without that step, property
uncertainty would propagate into nothing and the study would be theatre.

Bands are deliberately generous rather than optimistic: `mu ±20%`, `k ±15%`, `rho ±2%`,
`cp ±3%`.

| Property | Band | Swing in return temperature |
|---|---|---|
| `k` | ±15 % | **15.28 K** |
| `mu` | ±20 % | **14.37 K** |
| `cp` | ±3 % | 0.21 K |
| `rho` | ±2 % | 0.00 K |

All 17 corner combinations solve. `T_return` spans **850.67 – 880.69 K** (spread 30.0 K)
about a nominal 864.39 K, and the minimum freeze margin across every corner is 118.7 K.

**The finding worth carrying forward:** conductivity and viscosity are where the model's
real uncertainty lives; density and heat capacity are irrelevant at these bands. The 30 K
spread is the number to quote alongside any FLiBe result, not the 864.39 K point value.

### 3.1 How much of this is extrapolation? Less than §2 implies

Checked rather than assumed. The nominal loop spans **864.4 K (return) to 918.1 K (core
outlet)**:

- `k` (fitted 800–1200 K): **fully inside** its fitted range. Not extrapolated at all.
- `mu` (fitted 873–1073 K): inside at the hot end; the cold leg dips **8.6 K** below the
  873 K floor. At the property-corner extremes (850.7–934.4 K) that grows to ~22 K.

So at this operating point the correlations are barely extrapolated, and the 30 K spread
comes from the ±15%/±20% literature-scatter bands imposed on them — not from extrapolation
error. That is a materially different situation from the companion stellarator CFD study,
which pushes the same correlations across the whole 732–1703 K liquid range.

This corrects an overstatement earlier in this document's own §2: the extrapolation concern
is real for the correlations in general, but small for *this* model at *this* operating
point.

### 3.2 What is robust vs what is not

| Quantity | Depends on | Uncertainty |
|---|---|---|
| core temperature rise, 53.75 K | `cp` only (`Q/(mdot*cp)`) | **±1.61 K** — the most robust number in the FLiBe model |
| return temperature, HX sizing, wall temperatures | `k`, `mu` through the heat-transfer coefficient | **±15 K** |

The design margin survives either way: no corner comes within 118 K of freezing.

---

## 4. Honest summary

The equations are solved correctly — established three ways. The steam-side equations are
the right ones. The FLiBe-side equations are defensible but rest on extrapolated
correlations, and the resulting uncertainty is quantified rather than assumed away. Nothing
here has been compared against a physical experiment, so no result in this repository should
be presented as validated against measurement.

---

## Appendix — falsifiability rules referenced by the tests

| Rule | Requirement |
|---|---|
| **R1** Deletion test | Every physics term has a named test that fails if the term is deleted (`FALSIFIABILITY.md`). |
| **R2** No definitional balances | Every quantity in a conservation check is computed independently, e.g. `Q_hx` comes from the heat-exchanger model, never from `Q_core - W_pump`. |
| **R3** No dead exports | Every public symbol is exercised by a test that fails if its body is replaced with `raise NotImplementedError`. |
| **R4** Convergence is proven | Iterative solvers compute an explicit residual and raise if it is not met, rather than returning the last iterate. |
| **R5** Sensitivity is mandatory | Each governing input (`UA`, `mdot`, `tau`, `Q`) has a test asserting the direction and size of its effect. |
