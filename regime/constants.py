"""Physical constants and air species data."""

R_UNIV = 8.314462618      # J/(mol K)
K_B = 1.380649e-23         # J/K
N_A = 6.02214076e23        # 1/mol
AMU = 1.66053907e-27       # kg
ATM = 101325.0             # Pa

# Diatomic air (Ar lumped into N2). theta_v [K], M [kg/mol], mole fraction
SPECIES = {
    "N2": {"theta_v": 3393.0, "M": 28.0134e-3, "X": 0.79},
    "O2": {"theta_v": 2273.0, "M": 31.9988e-3, "X": 0.21},
}

# Mixture molar mass / gas constant (frozen composition)
M_MIX = sum(s["M"] * s["X"] for s in SPECIES.values())
R_MIX = R_UNIV / M_MIX

# Mass fractions
for _s in SPECIES.values():
    _s["Y"] = _s["M"] * _s["X"] / M_MIX
    _s["R"] = R_UNIV / _s["M"]
