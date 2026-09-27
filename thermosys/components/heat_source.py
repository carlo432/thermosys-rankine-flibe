def heat_source(inlet, mdot, Q_total):
    f=inlet.fluid; T=inlet.T+Q_total/(mdot*f.cp(inlet.T,inlet.P)); return f.state_TP(T,inlet.P)
