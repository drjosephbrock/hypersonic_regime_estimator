"""Text output: single-point report and sweep table."""


def _fmt_da(d):
    return "  ".join(f"{k}={v:.2e}" for k, v in d.items()) or "-"


def single(r):
    lines = [
        f"Freestream   alt={r.alt_km:.1f} km  V={r.V:.0f} m/s  M={r.mach:.2f}",
        f"             T1={r.T1:.1f} K  p1={r.p1:.4g} Pa  rho1={r.rho1:.4g} kg/m3",
        f"Total enth.  h0={r.h0/1e6:.3f} MJ/kg",
        f"Post-shock   T2_frozen={r.T2_frozen:.0f} K  T2_eq(vib)={r.T2_eq:.0f} K  "
        f"p2={r.p2:.4g} Pa  u2={r.u2:.0f} m/s",
        f"             e_vib/h at T2_eq = {100*r.vib_frac:.1f}%",
        f"Flow time    L={r.L:.3e} m  t_flow={r.t_flow:.3e} s",
        f"Da_vib       {_fmt_da(r.da_vib)}   -> {r.vib_regime}",
        f"Da_chem      {_fmt_da(r.da_chem)}   -> {r.chem_regime}",
        "",
        f"RECOMMEND:   {r.model}",
    ]
    lines += [f"  note: {n}" for n in r.notes]
    return "\n".join(lines)


def table(results):
    hdr = (f"{'alt km':>7} {'V m/s':>7} {'M':>6} {'h0 MJ/kg':>9} {'T2eq K':>7} "
           f"{'Da_vN2':>9} {'Da_cO2':>9} {'vib':>14} {'chem':>14}  model")
    out = [hdr, "-" * len(hdr)]
    for r in results:
        dc = r.da_chem.get("O2")
        out.append(
            f"{r.alt_km:7.1f} {r.V:7.0f} {r.mach:6.2f} {r.h0/1e6:9.3f} {r.T2_eq:7.0f} "
            f"{r.da_vib['N2'] if r.da_vib else float('nan'):9.2e} "
            f"{(f'{dc:9.2e}' if dc is not None else '        -')} "
            f"{r.vib_regime:>14} {r.chem_regime:>14}  {r.model}"
        )
    return "\n".join(out)
