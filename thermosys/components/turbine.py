import warnings
def turbine(inlet, P_out, eta_s):
    hs = inlet.fluid.state_PS(P_out, inlet.s).h
    h = inlet.h - eta_s*(inlet.h-hs)
    out = inlet.fluid.state_PH(P_out,h)
    if out.x is not None and out.x < 0.88: warnings.warn(f"turbine exit quality {out.x:.6f} < 0.88", UserWarning, stacklevel=2)
    return out, inlet.h-h
