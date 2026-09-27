"""Flow exergy and component exergy destruction.

Dead state 298.15 K / 101325 Pa.  psi = (h - h0) - T0*(s - s0); destruction = T0*s_gen.
"""
T0_DEAD = 298.15
P0_DEAD = 101325.0


def dead_state(fluid, T0=T0_DEAD, P0=P0_DEAD):
    return fluid.state_TP(T0, P0)


def flow_exergy(state, T0=T0_DEAD, P0=P0_DEAD):
    ref = dead_state(state.fluid, T0, P0)
    return (state.h - ref.h) - T0 * (state.s - ref.s)


def destruction(inlet, outlet, T0=T0_DEAD):
    """Adiabatic component: X_dest = T0 * (s_out - s_in), per unit mass."""
    return T0 * (outlet.s - inlet.s)


def turbine_second_law_efficiency(inlet, outlet, w_actual, T0=T0_DEAD, P0=P0_DEAD):
    drop = flow_exergy(inlet, T0, P0) - flow_exergy(outlet, T0, P0)
    return w_actual / drop


def pump_second_law_efficiency(inlet, outlet, w_actual, T0=T0_DEAD, P0=P0_DEAD):
    rise = flow_exergy(outlet, T0, P0) - flow_exergy(inlet, T0, P0)
    return rise / w_actual


def cycle_exergy_report(result, T0=T0_DEAD, P0=P0_DEAD):
    """Exergy accounting for a simple Rankine result dict (states 1-4)."""
    s1, s2, s3, s4 = result["states"][:4]
    return dict(
        psi=[flow_exergy(s, T0, P0) for s in (s1, s2, s3, s4)],
        x_dest_pump=destruction(s1, s2, T0),
        x_dest_turbine=destruction(s3, s4, T0),
        eta_II_turbine=turbine_second_law_efficiency(s3, s4, result["w_turbine"], T0, P0),
        eta_II_pump=pump_second_law_efficiency(s1, s2, result["w_pump"], T0, P0),
    )
