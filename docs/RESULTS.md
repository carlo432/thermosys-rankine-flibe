# Detailed results

All values regenerate from the scripts in `runs/` and are asserted by the test suite.

### Cycles

| Cycle | eta_th | exit quality |
|---|---:|---:|
| A: ideal 3 MPa → 75 kPa | 0.248878 | 0.796538 |
| B: 8 MPa/500 °C → 10 kPa, ηp=.85 ηt=.88 | 0.346279 | 0.874033 |
| C: ideal counterpart | 0.394173 | 0.810359 |
| Reheat, ideal | 0.408355 | 0.948708 |
| Reheat, ηp=.85 ηt=.90 | 0.372481 | 0.991263 |
| Regenerative, one open FWH (y = 0.214965) | 0.421198 | 0.810359 |

Reheat pressure optimum, swept: **1.844 MPa at η = 0.410002** — interior, as it should be.

### FLiBe loop (500 MW reference configuration)

`T_return = 903.691595 K`, `T_core_out = 957.445088 K` (rise 53.753493 K),
`W_pump = 197 754.677 W`, `Q_hx = 500 197 754.677 W`, `T_sec_out = 838.476750 K`,
`eps = 0.341411051`, `NTU = 0.483590173`, energy-balance residual 7.2e-16, freeze margin
171.69 K. UA sensitivity reproduces 959.211722 / 903.691595 / 875.960871 K.

### Coupled plant (UA = 2.5 MW/K)

| Quantity | Value |
|---|---:|
| `T_return` / `T_core_out` | 782.0236 / 835.7765 K |
| interface approach | 62.626 K |
| `mdot_steam` | 130.4050 kg/s |
| `Q_to_steam` | 500.192 MW |
| `W_net` | 204.064 MW |
| `eta_plant` vs standalone reheat | 0.408128 vs 0.408355 |
| penalty | 0.0227 efficiency points |
| freeze margin | 50.02 K |

### Transient scenarios (`runs/`, CSV + PNG each)

| Scenario | Result |
|---|---|
| +10% core power | settles at 914.0568 K, equals the algebraic solve at 550 MW |
| −20% mdot (dP = 0) | core rise ratio 1.250000 vs exact 3900/3120 = 1.250000 |
| pump coastdown | max core outlet 1093.85 K at t = 20 s |
| loss of heat sink | slope 1.496827 K/s vs analytic Q/(M·cp) = 1.496827 K/s |

### Moving-boundary boiler


Distributed enthalpy field (N finite volumes with real energy balances) with the phase
boundaries *diagnosed* by locating the h_f and h_g crossings and interpolating on the
downstream-face grid — which makes the diagnosis **exact and grid-independent** for the
linear steady profile (identical zone lengths at N = 10, 20, 60, 200). Deliberately not the
classical three-state moving-boundary ODE, which goes singular when a zone vanishes; here a
vanished zone reports zero length, and a test exercises exactly that. Bounding
simplification, stated not buried: uniform constant cell inventory. Superheat onset at
89.4 MW for the reference tube; transport delay matches the tube fill time M/mdot = 10 s.


## Modeling notes

1. **`ConstantCpFluid` added.** The loop specification describes the secondary stream by `(mdot, cp, T_in)`, but
   `counterflow_hx` takes `State` objects. This is the minimal `Fluid` that expresses that
   stream. It is not a third working fluid and carries no EOS.

2. **Plant energy-balance sign.** An early formulation wrote
   `Q_core = Q_to_steam + W_pump_flibe`. Pump work is added *to* the coolant, so the correct
   statement is `Q_to_steam = Q_core + W_pump`, which is consistent with the loop anchor values.

3. **The plant charges the primary for the reheater, not just the boiler.** Using
   `q_boiler` gave `eta_plant = 0.489 > eta_steam_alone = 0.408` — the cycle was being handed
   its reheat energy for free. The primary now supplies `q_in = q_boiler + q_reheat`, and the
   penalty comes out positive (0.023 points), as it must. Default `UA` is 2.5e6 W/K; the
   feasible band is roughly 1.5e6–4.3e6 W/K, outside which the return temperature leaves the
   liquid range and `LoopDomainError` is raised.

4. **Three test tolerances, with measurements:**
   - Isentropic exergy destruction: target was `< 1e-9` relative. The turbine path is
     exactly 0; the pump path floors at **6.6e-8** because CoolProp's `state_PS → state_PH`
     round trip carries ~2.7e-9 relative entropy error over compressed liquid. Library
     precision, not physics. Used 1e-6 for the pump, 1e-12 for the turbine.
   - FLiBe `mu` at 1703 K: the reference table prints 5 decimals there, so only 3 significant
     figures. `rtol=5e-3` for that row, 1e-3 elsewhere. 
   - Erlang response at `t = 1.0 s`: reference quotes `1.721e-04` to four significant figures;
     exact value is `1.7211562995589347e-04`. `rtol=1e-3` for that row, 1e-6 for the rest.

5. **`FlibeLoop.residual` evaluates bare temperatures with clamped properties** rather than
   building `State` objects. While bracketing, trial temperatures legitimately fall outside
   the liquid range; constructing states there raises and the scan then cannot tell which
   side of the range the root lies on — and so could not report the direction in
   `LoopDomainError`. `counterflow_hx` still produces the final states, and
   the mutation `loop duty from the exchanger` confirms the duty is genuinely the
   exchanger's.

6. **Coastdown "peak" is the maximum within the 20 s window**, not a turning point — the core
   outlet is still rising when the window closes. Labelled that way in the figure and the
   CSV rather than implying a peak was found.

7. **`FlibeLoopTransient.cp`** is evaluated at the coolant's own `T_min`, not at the
   secondary inlet temperature, which can legitimately sit outside the FLiBe liquid range.

