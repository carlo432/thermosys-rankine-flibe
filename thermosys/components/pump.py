def pump(inlet, P_out, eta_s):
    hs = inlet.fluid.state_PS(P_out, inlet.s).h
    h = inlet.h + (hs-inlet.h)/eta_s
    out = inlet.fluid.state_PH(P_out,h)
    return out, h-inlet.h
