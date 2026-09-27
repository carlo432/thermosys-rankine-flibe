"""Reheat Rankine cycle.

pump -> boiler -> HP turbine -> reheater -> LP turbine -> condenser.
Reheat exists to keep the LP exhaust dry; `x6` is therefore a headline result, not a note.
"""
from thermosys.components import pump, turbine, boiler, condenser
from thermosys.fluids import WaterFluid


class ReheatCycle:
    def __init__(self, fluid=None, P_boiler=8e6, P_reheat=1e6, P_cond=1e4,
                 T_boiler=773.15, T_reheat=773.15, eta_p=1.0, eta_t=1.0):
        self.fluid = fluid or WaterFluid()
        self.P_boiler, self.P_reheat, self.P_cond = P_boiler, P_reheat, P_cond
        self.T_boiler, self.T_reheat = T_boiler, T_reheat
        self.eta_p, self.eta_t = eta_p, eta_t

    def solve(self):
        f = self.fluid
        s1 = f.state_PQ(self.P_cond, 0.0)
        s2, w_p = pump(s1, self.P_boiler, self.eta_p)
        s3, q_boiler = boiler(s2, self.P_boiler, self.T_boiler, superheat=True)
        s4, w_hp = turbine(s3, self.P_reheat, self.eta_t)
        s5 = f.state_TP(self.T_reheat, self.P_reheat)
        q_reheat = s5.h - s4.h
        s6, w_lp = turbine(s5, self.P_cond, self.eta_t)

        q_in = q_boiler + q_reheat
        q_out = s6.h - s1.h
        w_t = w_hp + w_lp
        w_net = w_t - w_p
        return dict(
            states=(s1, s2, s3, s4, s5, s6),
            h1=s1.h, h2=s2.h, h3=s3.h, h4=s4.h, h5=s5.h, h6=s6.h,
            q_boiler=q_boiler, q_reheat=q_reheat, q_in=q_in, q_out=q_out,
            w_turbine=w_t, w_hp=w_hp, w_lp=w_lp, w_pump=w_p, w_net=w_net,
            eta_th=w_net / q_in, x6=s6.x, back_work_ratio=w_p / w_t,
        )
