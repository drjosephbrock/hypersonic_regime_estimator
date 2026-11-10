"""Combine freestream, shock, and relaxation times into a regime call."""

from dataclasses import dataclass, field

from . import thermo, relaxation as rx
from .atmosphere import us76
from .constants import R_MIX
from .shock import normal_shock

# --- thresholds (edit to taste) ---
T_VIB_ONSET = 800.0      # K, vibration starts to matter
T_O2_DISS = 2500.0       # K, O2 dissociation onset (~1 atm; lower at low p)
T_N2_DISS = 4000.0       # K, N2 dissociation onset
DA_EQ = 10.0             # Da > this -> equilibrium
DA_FROZEN = 0.1          # Da < this -> frozen
VIB_FRAC_MIN = 0.05      # e_vib/h below this -> vib nonequilibrium energetically minor


@dataclass
class Result:
    alt_km: float
    V: float
    T1: float
    p1: float
    rho1: float
    mach: float
    h0: float                 # J/kg
    T2_frozen: float          # vib frozen behind shock
    T2_eq: float              # vib equilibrated (no chemistry)
    p2: float
    u2: float
    vib_frac: float           # e_v/h at T2_eq
    L: float                  # flow length used
    t_flow: float
    da_vib: dict = field(default_factory=dict)
    da_chem: dict = field(default_factory=dict)
    vib_regime: str = ""
    chem_regime: str = ""
    model: str = ""
    notes: list = field(default_factory=list)


def _bin(da):
    if da >= DA_EQ:
        return "equilibrium"
    if da <= DA_FROZEN:
        return "frozen"
    return "nonequilibrium"


def _bin_all(das):
    """Regime over several species: uniform -> that regime, mixed -> nonequilibrium."""
    bins = {_bin(d) for d in das}
    return bins.pop() if len(bins) == 1 else "nonequilibrium"


def standoff(R_nose, eps):
    """Sphere shock standoff, Billig-type: delta ~ 0.78 R eps."""
    return 0.78 * R_nose * eps


def evaluate(alt_km, V, R_nose=0.1, L=None, T1=None, p1=None):
    """Classify flow regime at stagnation region.

    alt_km   : altitude (US76) unless T1/p1 given
    V        : velocity [m/s]
    R_nose   : nose radius [m], sets standoff/flow length
    L        : override flow length [m]
    """
    if T1 is None or p1 is None:
        T1, p1, rho1 = us76(alt_km * 1e3)
    else:
        rho1 = p1 / (R_MIX * T1)

    mach = V / (thermo.gamma(T1) * R_MIX * T1) ** 0.5
    h0 = thermo.enthalpy(T1) + 0.5 * V * V
    fr = normal_shock(T1, p1, V, vib=False)
    eq = normal_shock(T1, p1, V, vib=True)

    if L is None:
        L = standoff(R_nose, eq.eps)
    t_flow = rx.flow_time(L, eq.u)

    r = Result(alt_km, V, T1, p1, rho1, mach, h0, fr.T, eq.T, eq.p, eq.u,
               thermo.vib_fraction(eq.T), L, t_flow)

    # Vibration: evaluate tau at both ends of relaxation, keep the slower (conservative)
    for sp in ("N2", "O2"):
        tau = max(rx.tau_vib(sp, fr.T, fr.p), rx.tau_vib(sp, eq.T, eq.p))
        r.da_vib[sp] = rx.damkohler(t_flow, tau)

    # Chemistry: evaluate at T2_eq (upper bound w/o chemistry); report if hot enough
    if eq.T > T_O2_DISS:
        r.da_chem["O2"] = rx.damkohler(t_flow, rx.tau_diss("O2", eq.T, eq.p))
    if eq.T > T_N2_DISS:
        r.da_chem["N2"] = rx.damkohler(t_flow, rx.tau_diss("N2", eq.T, eq.p))

    _decide(r)
    return r


def _decide(r):
    if r.T2_eq < T_VIB_ONSET:
        r.vib_regime = "negligible"
        r.chem_regime = "none"
        r.model = "Calorically perfect gas (gamma = 1.4)"
        return

    r.vib_regime = _bin_all(r.da_vib.values())
    r.chem_regime = _bin_all(r.da_chem.values()) if r.da_chem else "none"

    v, c = r.vib_regime, r.chem_regime
    if c == "frozen":  # reactions too slow to matter in the shock layer
        c = "none"
    if c == "none" and r.vib_frac < VIB_FRAC_MIN:
        r.model = (f"Perfect gas / 1-T thermally perfect adequate "
                   f"(e_vib/h < {100*VIB_FRAC_MIN:.0f}%); 2-T only for fidelity")
        return
    if v == "equilibrium":
        if c == "none":
            r.model = "1-T thermally perfect gas (vib in equilibrium, no chemistry)"
        elif c == "equilibrium":
            r.model = "1-T equilibrium air (equilibrium chemistry)"
        else:
            r.model = "1-T finite-rate chemistry (vib in equilibrium)"
    elif v == "nonequilibrium":
        r.model = "2-T (T-Tv) model" + (", finite-rate chemistry" if c != "none" else "")
    elif c == "none":  # vib frozen, no chemistry
        r.model = ("Perfect gas OK in shock layer (vib frozen); "
                   "2-T if downstream relaxation (long body/wake) matters")
    else:  # vib frozen, chemistry active
        r.model = "2-T model, finite-rate chemistry (vib frozen in shock layer)"

    if r.T2_eq > T_O2_DISS:
        r.notes.append("T2_eq ignores dissociation -> overestimates T; Da_chem is an upper bound.")
    if r.p2 < 1e3 and r.T2_eq > T_O2_DISS:
        r.notes.append("Low p2: dissociation onsets shift lower than the ~1 atm thresholds.")
    if r.alt_km > 70:
        r.notes.append("Check Knudsen number; continuum may be marginal.")
