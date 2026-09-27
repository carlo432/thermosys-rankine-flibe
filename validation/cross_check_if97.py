#!/usr/bin/env python3
"""Independent cross-check of the cycle anchors against IAPWS-IF97.

Deliberately imports NOTHING from `thermosys`. The property backend is different
(IF97 industrial formulation vs CoolProp's IF95 Helmholtz) and the cycle algebra is
written a second time from scratch, so agreement is evidence about the cycles rather
than a restatement of them.

Scope limit, stated up front: IF97 and IF95 agree to roughly 0.01-0.1% by construction.
This check can catch a wrong state point, a flipped isentropic relation, or a mis-assembled
efficiency. It CANNOT confirm the anchors at the 1e-6 tolerance the test suite asserts.

Run:  python3 validation/cross_check_if97.py
"""
from iapws import IAPWS97

# IF97 works in MPa / kJ per kg; the anchors are SI. Convert at the boundary only.
MPa = 1e-6
kJ = 1e3


def sat(P, x):
    return IAPWS97(P=P * MPa, x=x)


def at_PT(P, T):
    return IAPWS97(P=P * MPa, T=T)


def at_Ps(P, s):
    return IAPWS97(P=P * MPa, s=s)


def at_Ph(P, h):
    return IAPWS97(P=P * MPa, h=h)


def pump(state_in, P_out, eta):
    ideal = at_Ps(P_out, state_in.s)
    h = state_in.h + (ideal.h - state_in.h) / eta
    return at_Ph(P_out, h), (h - state_in.h) * kJ


def turbine(state_in, P_out, eta):
    ideal = at_Ps(P_out, state_in.s)
    h = state_in.h - eta * (state_in.h - ideal.h)
    return at_Ph(P_out, h), (state_in.h - h) * kJ


def rankine(P_boiler, P_cond, T_boiler=None, eta_p=1.0, eta_t=1.0):
    s1 = sat(P_cond, 0.0)
    s2, w_p = pump(s1, P_boiler, eta_p)
    s3 = at_PT(P_boiler, T_boiler) if T_boiler else sat(P_boiler, 1.0)
    s4, w_t = turbine(s3, P_cond, eta_t)
    q_in = (s3.h - s2.h) * kJ
    return dict(eta_th=(w_t - w_p) / q_in, x4=s4.x, q_in=q_in, w_t=w_t, w_p=w_p)


def reheat(P_boiler, P_reheat, P_cond, T_hot, eta_p=1.0, eta_t=1.0):
    s1 = sat(P_cond, 0.0)
    s2, w_p = pump(s1, P_boiler, eta_p)
    s3 = at_PT(P_boiler, T_hot)
    s4, w_hp = turbine(s3, P_reheat, eta_t)
    s5 = at_PT(P_reheat, T_hot)
    s6, w_lp = turbine(s5, P_cond, eta_t)
    q_in = ((s3.h - s2.h) + (s5.h - s4.h)) * kJ
    return dict(eta_th=(w_hp + w_lp - w_p) / q_in, x6=s6.x, q_in=q_in,
                w_t=w_hp + w_lp, w_p=w_p)


def regenerative(P_boiler, P_extract, P_cond, T_hot):
    s1 = sat(P_cond, 0.0)
    s2, w_p1 = pump(s1, P_extract, 1.0)
    s3 = sat(P_extract, 0.0)
    s4, w_p2 = pump(s3, P_boiler, 1.0)
    s5 = at_PT(P_boiler, T_hot)
    s6, w_hp = turbine(s5, P_extract, 1.0)
    s7, w_lp = turbine(s6, P_cond, 1.0)
    y = (s3.h - s2.h) / (s6.h - s2.h)
    q_in = (s5.h - s4.h) * kJ
    w_t = w_hp + (1.0 - y) * w_lp
    w_p = (1.0 - y) * w_p1 + w_p2
    return dict(eta_th=(w_t - w_p) / q_in, y=y, q_in=q_in, w_t=w_t, w_p=w_p)


CASES = [
    ("A  ideal 3 MPa sat -> 75 kPa", lambda: rankine(3e6, 75e3), "eta_th", 0.248878),
    ("B  8 MPa/500C -> 10 kPa, .85/.88",
     lambda: rankine(8e6, 1e4, 773.15, 0.85, 0.88), "eta_th", 0.346279),
    ("C  ideal 8 MPa/500C -> 10 kPa",
     lambda: rankine(8e6, 1e4, 773.15), "eta_th", 0.394173),
    ("R  reheat ideal, 1 MPa reheat",
     lambda: reheat(8e6, 1e6, 1e4, 773.15), "eta_th", 0.408355),
    ("R  reheat .85/.90",
     lambda: reheat(8e6, 1e6, 1e4, 773.15, 0.85, 0.90), "eta_th", 0.372481),
    ("G  regenerative, 1 MPa extraction",
     lambda: regenerative(8e6, 1e6, 1e4, 773.15), "eta_th", 0.421198),
    ("G  extraction fraction y",
     lambda: regenerative(8e6, 1e6, 1e4, 773.15), "y", 0.214965),
    ("A  turbine exit quality x4", lambda: rankine(3e6, 75e3), "x4", 0.796538),
    ("C  turbine exit quality x4",
     lambda: rankine(8e6, 1e4, 773.15), "x4", 0.810359),
    ("R  LP exit quality x6", lambda: reheat(8e6, 1e6, 1e4, 773.15), "x6", 0.948708),
]

SAT_ANCHORS = [
    (1e4, 318.9563, 191805.94, 2583858.67),
    (75e3, 364.9080, 384437.44, 2662386.22),
    (101325.0, 373.1243, 419057.73, 2675529.33),
    (3e6, 507.0031, 1008344.61, 2803153.10),
    (8e6, 568.1577, 1317311.05, 2758683.20),
]


def main():
    print("IAPWS-IF97 cross-check  (independent formulation and independent cycle code)")
    print("=" * 78)
    print("\nSaturation anchors")
    print(f"  {'P [Pa]':>10} {'T_sat dev':>12} {'h_f dev':>12} {'h_g dev':>12}")
    worst_sat = 0.0
    for P, T_ref, hf_ref, hg_ref in SAT_ANCHORS:
        f, g = sat(P, 0.0), sat(P, 1.0)
        dT = (f.T - T_ref) / T_ref
        dhf = (f.h * kJ - hf_ref) / hf_ref
        dhg = (g.h * kJ - hg_ref) / hg_ref
        worst_sat = max(worst_sat, abs(dT), abs(dhf), abs(dhg))
        print(f"  {P:10.0f} {dT:11.2e} {dhf:11.2e} {dhg:11.2e}")

    print("\nCycle anchors")
    print(f"  {'case':36} {'IF97':>10} {'anchor':>10} {'deviation':>11}")
    worst = 0.0
    for label, fn, key, ref in CASES:
        got = fn()[key]
        dev = (got - ref) / ref
        worst = max(worst, abs(dev))
        flag = "" if abs(dev) < 2e-3 else "   <-- CHECK"
        print(f"  {label:36} {got:10.6f} {ref:10.6f} {dev:10.2e}{flag}")

    print("\n" + "=" * 78)
    print(f"worst saturation deviation : {worst_sat:.2e}")
    print(f"worst cycle deviation      : {worst:.2e}")
    print("\nIF97 vs IF95 mutual tolerance is ~1e-4..1e-3 on enthalpy, so deviations at")
    print("that level are the two formulations differing, not an error in either. A wrong")
    print("state point or a flipped isentropic relation would show up orders of magnitude")
    print("larger. This is verification of the cycle construction, not validation.")
    return 0 if worst < 5e-3 else 1


if __name__ == "__main__":
    raise SystemExit(main())
