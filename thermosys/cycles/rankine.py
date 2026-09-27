from thermosys.components import pump,turbine,boiler,condenser
class RankineCycle:
    def __init__(self,fluid=None,P_boiler=3e6,P_cond=75e3,T_boiler=None,eta_p=1.,eta_t=1.):
        from thermosys.fluids import WaterFluid
        self.fluid=fluid or WaterFluid(); self.P_boiler=P_boiler; self.P_cond=P_cond; self.T_boiler=T_boiler; self.eta_p=eta_p; self.eta_t=eta_t
    def solve(self):
        f=self.fluid; s1=f.state_PQ(self.P_cond,0); s2,wp=pump(s1,self.P_boiler,self.eta_p)
        s3,q=boiler(s2,self.P_boiler,self.T_boiler,self.T_boiler is not None)
        s4,wt=turbine(s3,self.P_cond,self.eta_t)
        qout=s4.h-s1.h; wnet=wt-wp
        return dict(states=(s1,s2,s3,s4),h1=s1.h,h2=s2.h,h3=s3.h,h4=s4.h,h2s=f.state_PS(self.P_boiler,s1.s).h,h4s=f.state_PS(self.P_cond,s3.s).h,q_in=q,q_out=qout,w_pump=wp,w_turbine=wt,w_net=wnet,eta_th=wnet/q,back_work_ratio=wp/wt)

def rankine_cycle(**kwargs):
    return RankineCycle(**kwargs).solve()
