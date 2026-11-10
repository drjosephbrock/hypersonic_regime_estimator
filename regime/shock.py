"""Normal-shock jump for thermally perfect air (no chemistry)."""

from dataclasses import dataclass

from . import thermo
from .constants import R_MIX


@dataclass
class ShockState:
    T: float      # K
    p: float      # Pa
    rho: float    # kg/m^3
    u: float      # m/s
    eps: float    # rho1/rho2


def normal_shock(T1, p1, u1, vib=True, tol=1e-10, max_iter=200):
    """Jump across normal shock.

    vib=False -> vibration frozen at T1 (state immediately behind shock).
    vib=True  -> vibration equilibrated (thermally perfect, no dissociation).
    """
    rho1 = p1 / (R_MIX * T1)
    h1 = thermo.enthalpy(T1, vib)
    eps, T2 = 0.2, T1
    for _ in range(max_iter):
        p2 = p1 + rho1 * u1 * u1 * (1.0 - eps)
        h2 = h1 + 0.5 * u1 * u1 * (1.0 - eps * eps)
        T2 = thermo.temperature(h2, vib, T0=T2)
        rho2 = p2 / (R_MIX * T2)
        eps_new = rho1 / rho2
        if abs(eps_new - eps) < tol:
            eps = eps_new
            break
        eps = 0.5 * (eps + eps_new)  # relax for robustness
    else:
        raise RuntimeError("normal shock did not converge")
    return ShockState(T2, p2, rho2, u1 * eps, eps)
