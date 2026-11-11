"""Optional V-altitude regime map (requires matplotlib + numpy)."""


def regime_map(alts_km, vels, R_nose=0.1, fname="regime_map.png"):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    from .classify import evaluate

    # one class per (vib, chem) model from classify._decide, so chemistry is never hidden
    labels = ["Perfect gas", "1-T thermally perfect", "2-T, no chemistry",
              "2-T + finite-rate chem", "1-T + finite-rate chem", "1-T equilibrium air"]
    Z = np.zeros((len(alts_km), len(vels)))
    H0 = np.zeros_like(Z)
    for i, a in enumerate(alts_km):
        print(i,len(alts_km),a)
        for j, v in enumerate(vels):
            r = evaluate(a, v, R_nose=R_nose)
            H0[i, j] = r.h0 / 1e6
            m = r.model
            if m.startswith(("Calorically", "Perfect")):
                Z[i, j] = 0
            elif m.startswith("1-T thermally"):
                Z[i, j] = 1
            elif m.startswith("2-T"):
                Z[i, j] = 3 if "chemistry" in m else 2
            elif m.startswith("1-T finite-rate"):
                Z[i, j] = 4
            else:  # 1-T equilibrium air
                Z[i, j] = 5

    cmap = ListedColormap(["#d9e6f2", "#9ec5e8", "#e07a5f", "#9b2d20", "#f2c879", "#7fb069"])
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pcolormesh(np.asarray(vels) / 1e3, alts_km, Z, cmap=cmap, vmin=-0.5, vmax=5.5,
                  shading="nearest")
    cs = ax.contour(np.asarray(vels) / 1e3, alts_km, H0, levels=[1, 2.5, 5, 10, 20],
                    colors="k", linewidths=0.7)
    ax.clabel(cs, fmt="h0=%g MJ/kg", fontsize=8)
    ax.set_xlabel("V [km/s]")
    ax.set_ylabel("Altitude [km]")
    ax.set_title(f"Thermochemical regime, stagnation region, R_nose={R_nose} m")
    import matplotlib.patches as mp
    present = set(np.unique(Z).astype(int))
    ax.legend(handles=[mp.Patch(color=cmap.colors[k], label=labels[k])
                       for k in range(len(labels)) if k in present],
              loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    return fname
