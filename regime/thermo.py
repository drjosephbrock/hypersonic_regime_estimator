"""Thermally perfect air (frozen composition, harmonic-oscillator vibration)."""

import math

from .constants import SPECIES, R_MIX


def e_vib_species(T, theta, R):
    """Vibrational energy per unit mass of one species [J/kg]."""
    x = theta / T
    return R * theta / math.expm1(x) if x < 700 else 0.0


def e_vib(T):
    """Mixture vibrational energy [J/kg]."""
    return sum(s["Y"] * e_vib_species(T, s["theta_v"], s["R"]) for s in SPECIES.values())


def h_tr(T):
    """Translational + rotational enthalpy (incl. pv) [J/kg]."""
    return 3.5 * R_MIX * T


def enthalpy(T, vib=True):
    """Static enthalpy [J/kg], zero at 0 K. vib=False -> calorically perfect."""
    return h_tr(T) + (e_vib(T) if vib else 0.0)


def cp(T, vib=True):
    """Specific heat at constant pressure [J/(kg K)]."""
    c = 3.5 * R_MIX
    if vib:
        for s in SPECIES.values():
            x = s["theta_v"] / T
            if x < 700:
                ex = math.exp(x)
                c += s["Y"] * s["R"] * x * x * ex / (ex - 1.0) ** 2
    return c


def gamma(T, vib=True):
    c = cp(T, vib)
    return c / (c - R_MIX)


def temperature(h, vib=True, T0=1000.0):
    """Invert enthalpy for T via Newton."""
    if not vib:
        return h / (3.5 * R_MIX)
    T = max(T0, 50.0)
    for _ in range(100):
        dT = (enthalpy(T) - h) / cp(T)
        T -= dT
        T = max(T, 10.0)
        if abs(dT) < 1e-8 * T:
            return T
    raise RuntimeError("enthalpy inversion did not converge")


def vib_fraction(T):
    """e_vib / h at T."""
    return e_vib(T) / enthalpy(T)
