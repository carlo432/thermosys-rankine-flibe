import math, warnings
from scipy.optimize import brentq
from thermosys.state import State

class FlibeFluid:
    name = "FLiBe"; T_min = 732.0; T_max = 1703.0; T_ref = 732.0
    def _check(self,T):
        if not self.T_min <= T <= self.T_max: raise ValueError(f"FLiBe temperature {T} K outside 732–1703 K")
    def _warn(self,T,lo,hi,prop):
        if not lo <= T <= hi: warnings.warn(f"{prop} correlation extrapolated at {T} K", UserWarning, stacklevel=3)
    def rho(self,T,P=0): self._check(T); return 2415.6-0.49072*T
    def mu(self,T,P=0): self._check(T); self._warn(T,873,1073,"mu"); return 1.16e-4*math.exp(3755.0/T)
    def cp(self,T,P=0): self._check(T); return 2386.0
    def k(self,T,P=0): self._check(T); self._warn(T,800,1200,"k"); return 1.212+9e-5*T
    def Pr(self,T,P=0): return self.mu(T,P)*self.cp(T,P)/self.k(T,P)
    def _s(self,T): return self.cp(T,0)*math.log(T/self.T_ref)
    def state_TP(self,T,P):
        self._check(T); return State(self,T,P,self.cp(T,P)*(T-self.T_ref),self._s(T),None)
    def state_PH(self,P,h): return self.state_TP(self.T_ref+h/self.cp(self.T_ref,P),P)
    def state_PS(self,P,s): return self.state_TP(self.T_ref*math.exp(s/self.cp(self.T_ref,P)),P)
    def state_PQ(self,P,Q): raise NotImplementedError("FLiBe has no two-phase quality")
    def T_from_nu(self,nu):
        return brentq(lambda T:self.mu(T,0)/self.rho(T,0)-nu,self.T_min,self.T_max)

FLiBeFluid = FlibeFluid
