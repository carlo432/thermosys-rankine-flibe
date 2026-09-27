#!/usr/bin/env python3
"""Rule R1 deletion check, run programmatically.

For each physics term below, neutralise it in the source, run the suite, and record which
tests fail.  A term whose deletion breaks nothing is dead weight and the run reports it as
SURVIVED -- that is a verification failure, not a curiosity.

Writes FALSIFIABILITY.md.  Restores every file it touches, including on error.

Run:  python3 tools/mutation_check.py
"""
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

# (label, relative path, snippet to replace, replacement)
MUTATIONS = [
    # -- fluid properties ---------------------------------------------------
    ("FLiBe density temperature dependence", "thermosys/fluids/flibe.py",
     "return 2415.6-0.49072*T", "return 2415.6-0.0*T"),
    ("FLiBe Arrhenius viscosity", "thermosys/fluids/flibe.py",
     "return 1.16e-4*math.exp(3755.0/T)", "return 1.16e-4*math.exp(3755.0/950.5)"),
    ("FLiBe conductivity temperature dependence", "thermosys/fluids/flibe.py",
     "return 1.212+9e-5*T", "return 1.1"),
    ("FLiBe liquid-range bound check", "thermosys/fluids/flibe.py",
     "if not self.T_min <= T <= self.T_max: raise ValueError",
     "if False: raise ValueError"),
    ("water single-phase quality mapped to None", "thermosys/fluids/water.py",
     "return None if x < 0.0 or x > 1.0 else x", "return x"),

    # -- components ---------------------------------------------------------
    ("pump isentropic efficiency direction", "thermosys/components/pump.py",
     "h = inlet.h + (hs-inlet.h)/eta_s", "h = inlet.h + (hs-inlet.h)*eta_s"),
    ("turbine isentropic efficiency", "thermosys/components/turbine.py",
     "h = inlet.h - eta_s*(inlet.h-hs)", "h = inlet.h - (inlet.h-hs)"),
    ("turbine wetness warning", "thermosys/components/turbine.py",
     "if out.x is not None and out.x < 0.88:", "if False:"),
    ("HX capacity-ratio term in effectiveness", "thermosys/components/hx.py",
     "return (1.0 - e) / (1.0 - Cr * e)", "return 1.0 - e"),
    ("HX Cr -> 1 limiting branch", "thermosys/components/hx.py",
     "if abs(Cr - 1.0) < 1e-12:\n        return NTU / (1.0 + NTU)",
     "if False:\n        return NTU / (1.0 + NTU)"),
    ("HX duty driving temperature difference", "thermosys/components/hx.py",
     "return eps * C_min * (T_hot_in - T_cold_in)", "return eps * C_min * 50.0"),

    # -- loop closure -------------------------------------------------------
    ("loop pump work in the energy balance", "thermosys/cycles/flibe_loop.py",
     "T_hot = T_return + (self.Q_core + W_pump) / C_h",
     "T_hot = T_return + self.Q_core / C_h"),
    ("loop duty from the exchanger (no definitional balance)",
     "thermosys/cycles/flibe_loop.py",
     "q = hx_duty(T_hot, self.T_sec_in, C_h, C_c, self.UA)",
     "q = self.Q_core + W_pump"),
    ("loop convergence residual check", "thermosys/cycles/flibe_loop.py",
     "if abs(r) > tol:", "if False:"),

    # -- cycles -------------------------------------------------------------
    ("reheater heat input", "thermosys/cycles/reheat.py",
     "q_in = q_boiler + q_reheat", "q_in = q_boiler"),
    ("extraction fraction from the FWH balance", "thermosys/cycles/regenerative.py",
     "y = (s3.h - s2.h) / (s6.h - s2.h)", "y = 0.2"),
    ("plant charges the primary for reheat duty", "thermosys/cycles/plant.py",
     'return r["states"][1], r["states"][2], r["q_in"]',
     'return r["states"][1], r["states"][2], r["q_boiler"]'),

    # -- transient right-hand sides -----------------------------------------
    ("single-CV advection term", "thermosys/transient/lumped.py",
     "return [(self.mdot * self.cp * (self.T_in - y[0]) + self.Q) / (self.M * self.cp)]",
     "return [(self.Q) / (self.M * self.cp)]"),
    ("single-CV heat input", "thermosys/transient/lumped.py",
     "+ self.Q) / (self.M * self.cp)]", "+ 0.0) / (self.M * self.cp)]"),
    ("CV-chain upstream coupling", "thermosys/transient/lumped.py",
     "upstream = np.concatenate(([1.0], y[:-1]))",
     "upstream = np.ones_like(y)"),
    ("loop network closure (np.roll)", "thermosys/transient/lumped.py",
     "upstream = np.roll(y, 1)", "upstream = y.copy()"),
    ("network core heat source", "thermosys/transient/lumped.py",
     "dTdt[self.i_core] += self.Q_at(t) + self.pump_work(y[-1], mdot)",
     "dTdt[self.i_core] += self.pump_work(y[-1], mdot)"),
    ("network HX heat removal", "thermosys/transient/lumped.py",
     "dTdt[self.i_hx] -= q", "dTdt[self.i_hx] -= 0.0"),
    ("network thermal capacitance", "thermosys/transient/lumped.py",
     "return dTdt / (self.masses * cp)", "return dTdt / (self.masses.mean() * cp)"),
    ("freeze/boil excursion check", "thermosys/transient/lumped.py",
     "if lo < self.fluid.T_min:", "if False:"),

    # -- exergy -------------------------------------------------------------
    ("dead-state temperature in exergy destruction", "thermosys/analysis/exergy.py",
     "return T0 * (outlet.s - inlet.s)", "return (outlet.s - inlet.s)"),

    # -- moving-boundary boiler ---------------------------------------------
    ("boiler enthalpy advection", "thermosys/transient/boiler_zones.py",
     "return (mdot * (upstream - h) + q_cell) / self.cell_mass",
     "return (0.0 * (upstream - h) + q_cell) / self.cell_mass"),
    ("boiler heat input", "thermosys/transient/boiler_zones.py",
     "+ q_cell) / self.cell_mass", "+ 0.0) / self.cell_mass"),
    ("boiler cell inventory", "thermosys/transient/boiler_zones.py",
     "return self.M_total / self.n_cells", "return self.M_total"),
    ("boiler downstream-face interpolation grid", "thermosys/transient/boiler_zones.py",
     "z = np.arange(self.n_cells + 1) * self.dz",
     "z = (np.arange(self.n_cells + 1) + 0.5) * self.dz"),
    ("boiler zone clipping to the tube", "thermosys/transient/boiler_zones.py",
     "L2 = min(max(L2, 0.0), self.length - L1)", "L2 = max(L2, 0.0)"),

    # -- heat-transfer correlation ------------------------------------------
    ("Dittus-Boelter Reynolds exponent", "thermosys/components/htc.py",
     "Nu = 0.023 * Re ** 0.8 * Pr", "Nu = 0.023 * Re ** 0.0 * Pr"),
    ("Dittus-Boelter conductivity scaling", "thermosys/components/htc.py",
     "return Nu * fluid.k(T, P) / D", "return Nu * 1.1 / D"),
    ("perturbation applied to Prandtl number", "thermosys/fluids/perturbed.py",
     "return self.mu(T, P) * self.cp(T, P) / self.k(T, P)",
     "return self.base.Pr(T, P)"),
]


def run_suite():
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--no-header", "--tb=no", "-p", "no:cacheprovider"],
        cwd=ROOT, capture_output=True, text=True, timeout=1800)
    failed = [ln.removeprefix("FAILED ").split(" ")[0].strip()
              for ln in proc.stdout.splitlines() if ln.startswith("FAILED ")]
    errored = any(ln.startswith("ERROR ") for ln in proc.stdout.splitlines())
    return failed, errored, proc.stdout.strip().splitlines()[-1] if proc.stdout else ""


def main():
    baseline_failed, _, baseline_line = run_suite()
    if baseline_failed:
        sys.exit(f"baseline suite is not green ({len(baseline_failed)} failures); fix first")
    print(f"baseline: {baseline_line}")

    rows, survivors = [], []
    for label, rel, old, new in MUTATIONS:
        path = ROOT / rel
        original = path.read_text()
        if old not in original:
            rows.append((label, rel, "SNIPPET NOT FOUND", ""))
            survivors.append(label)
            print(f"  !! {label}: snippet not found in {rel}")
            continue
        try:
            path.write_text(original.replace(old, new, 1))
            failed, errored, last = run_suite()
        finally:
            path.write_text(original)

        if failed or errored:
            killer = failed[0] if failed else "collection error"
            rows.append((label, rel, killer, str(len(failed) if failed else "collect")))
            print(f"  killed by {killer}  <- {label}")
        else:
            rows.append((label, rel, "SURVIVED", "0"))
            survivors.append(label)
            print(f"  !! SURVIVED: {label}")

    with open(ROOT / "FALSIFIABILITY.md", "w") as fh:
        fh.write("# Falsifiability map (rule R1)\n\n")
        fh.write("Generated by `python3 tools/mutation_check.py`. Each row names a physics "
                 "term, the mutation applied to it, and the first test that fails as a "
                 "result. Regenerate after changing the models.\n\n")
        fh.write(f"Baseline: `{baseline_line}`\n\n")
        fh.write("| Physics term | Module | Killed by | Failing tests |\n")
        fh.write("|---|---|---|---|\n")
        for label, rel, killer, n in rows:
            mark = "**SURVIVED**" if killer in ("SURVIVED", "SNIPPET NOT FOUND") else f"`{killer}`"
            fh.write(f"| {label} | `{rel}` | {mark} | {n} |\n")
        fh.write(f"\n{len(rows) - len(survivors)}/{len(rows)} mutations killed.\n")
        if survivors:
            fh.write("\n## Survivors — unfalsified terms\n\n")
            for s in survivors:
                fh.write(f"- {s}\n")
        else:
            fh.write("\nNo survivors: every term above is covered by at least one test "
                     "that fails when the term is removed.\n")

    print(f"\n{len(rows) - len(survivors)}/{len(rows)} killed -> FALSIFIABILITY.md")
    return 1 if survivors else 0


if __name__ == "__main__":
    sys.exit(main())
