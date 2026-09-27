"""Regenerative Rankine cycle with one open feedwater heater.

States: 1 condenser sat liquid -> 2 pump-1 to P_extract -> 3 FWH outlet (sat liquid at
P_extract) -> 4 pump-2 to P_boiler -> 5 boiler outlet -> 6 extraction -> 7 condenser inlet.

The extraction fraction y is SOLVED from the FWH energy balance, never hardcoded (rule R2):

    y*h6 + (1 - y)*h2 = h3   ->   y = (h3 - h2) / (h6 - h2)
"""
from thermosys.components import pump, turbine, boiler
from thermosys.fluids import WaterFluid


class RegenerativeCycle:
    def __init__(self, fluid=None, P_boiler=8e6, P_extract=1e6, P_cond=1e4,
                 T_boiler=773.15, eta_p=1.0, eta_t=1.0):
        self.fluid = fluid or WaterFluid()
        self.P_boiler, self.P_extract, self.P_cond = P_boiler, P_extract, P_cond
        self.T_boiler = T_boiler
        self.eta_p, self.eta_t = eta_p, eta_t

    def solve(self):
        f = self.fluid
        s1 = f.state_PQ(self.P_cond, 0.0)
        s2, w_p1 = pump(s1, self.P_extract, self.eta_p)
        s3 = f.state_PQ(self.P_extract, 0.0)          # FWH outlet: saturated liquid
        s4, w_p2 = pump(s3, self.P_boiler, self.eta_p)
        s5, _ = boiler(s4, self.P_boiler, self.T_boiler, superheat=True)
        s6, w_hp = turbine(s5, self.P_extract, self.eta_t)
        s7, w_lp = turbine(s6, self.P_cond, self.eta_t)

        y = (s3.h - s2.h) / (s6.h - s2.h)             # FWH energy balance, solved
        q_in = s5.h - s4.h
        q_out = (1.0 - y) * (s7.h - s1.h)
        w_t = w_hp + (1.0 - y) * w_lp
        w_p = (1.0 - y) * w_p1 + w_p2
        w_net = w_t - w_p
        return dict(
            states=(s1, s2, s3, s4, s5, s6, s7),
            h1=s1.h, h2=s2.h, h3=s3.h, h4=s4.h, h5=s5.h, h6=s6.h, h7=s7.h,
            y=y, T_fwh_sat=s3.T, q_in=q_in, q_out=q_out,
            w_turbine=w_t, w_pump=w_p, w_net=w_net,
            eta_th=w_net / q_in, x7=s7.x,
            fwh_residual=abs(y * s6.h + (1.0 - y) * s2.h - s3.h) / s3.h,
        )
