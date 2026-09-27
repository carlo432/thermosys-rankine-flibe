def boiler(inlet, P_out, T_out=None, superheat=False):
    f=inlet.fluid
    if superheat and T_out is not None: out=f.state_TP(T_out,P_out)
    else: out=f.state_PQ(P_out,1.0)
    return out, out.h-inlet.h
