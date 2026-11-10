"""Optional V-altitude regime map (requires matplotlib + numpy)."""


def regime_map(alts_km, vels, R_nose=0.1, fname="regime_map.png"):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    from .classify import evaluate

    codes = {"perfect": 0, "1T-TP": 1, "1T-chem": 2, "2T": 3}
    Z = np.zeros((len(alts_km), len(vels)))
    H0 = np.zeros_like(Z)
    for i, a in enumerate(alts_km):
        for j, v in enumerate(vels):
            r = evaluate(a, v, R_nose=R_nose)
            H0[i, j] = r.h0 / 1e6
            if r.vib_regime == "negligible" or r.model.startswith("Perfect"):
                Z[i, j] = codes["perfect"]
            elif r.model.startswith("2-T"):
                Z[i, j] = codes["2T"]
            elif r.chem_regime == "none":
                Z[i, j] = codes["1T-TP"]
            else:
                Z[i, j] = codes["1T-chem"]

    cmap = ListedColormap(["#d9e6f2", "#9ec5e8", "#f2c879", "#e07a5f"])
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pcolormesh(np.asarray(vels) / 1e3, alts_km, Z, cmap=cmap, vmin=-0.5, vmax=3.5,
                  shading="nearest")
    cs = ax.contour(np.asarray(vels) / 1e3, alts_km, H0, levels=[1, 2.5, 5, 10, 20],
                    colors="k", linewidths=0.7)
    ax.clabel(cs, fmt="h0=%g MJ/kg", fontsize=8)
    ax.set_xlabel("V [km/s]")
    ax.set_ylabel("Altitude [km]")
    ax.set_title(f"Thermochemical regime, stagnation region, R_nose={R_nose} m")
    import matplotlib.patches as mp
    ax.legend(handles=[mp.Patch(color=c, label=l) for c, l in zip(
        cmap.colors, ["Perfect gas", "1-T thermally perfect", "1-T chemistry", "2-T"])],
        loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    return fname
