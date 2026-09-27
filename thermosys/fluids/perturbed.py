"""A fluid whose properties are scaled by fixed factors.

Uncertainty propagation with no changes to any component or cycle: because everything
downstream talks to the Fluid protocol, perturbing a correlation is a wrapper, not an
edit. This is the architecture rule (§1) paying for itself.
"""


class PerturbedFluid:
    """Wraps a fluid, multiplying selected properties by constant factors."""

    def __init__(self, base, mu=1.0, k=1.0, rho=1.0, cp=1.0):
        self.base = base
        self.factors = dict(mu=mu, k=k, rho=rho, cp=cp)
        self.name = f"{base.name}(perturbed)"
        self.T_min, self.T_max = base.T_min, base.T_max

    def rho(self, T, P=0.0):
        return self.base.rho(T, P) * self.factors["rho"]

    def cp(self, T, P=0.0):
        return self.base.cp(T, P) * self.factors["cp"]

    def mu(self, T, P=0.0):
        return self.base.mu(T, P) * self.factors["mu"]

    def k(self, T, P=0.0):
        return self.base.k(T, P) * self.factors["k"]

    def Pr(self, T, P=0.0):
        return self.mu(T, P) * self.cp(T, P) / self.k(T, P)

    def state_TP(self, T, P):
        return self.base.state_TP(T, P)

    def state_PH(self, P, h):
        return self.base.state_PH(P, h)

    def state_PS(self, P, s):
        return self.base.state_PS(P, s)

    def state_PQ(self, P, Q):
        return self.base.state_PQ(P, Q)
