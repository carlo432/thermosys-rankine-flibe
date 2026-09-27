def condenser(inlet, P_out, subcool_K=0.0):
    f=inlet.fluid; sat=f.state_PQ(P_out,0.0)
    out = f.state_TP(sat.T-subcool_K,P_out) if subcool_K else sat
    return out, inlet.h-out.h
