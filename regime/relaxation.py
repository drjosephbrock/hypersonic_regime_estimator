"""Characteristic times: vibrational relaxation and dissociation.

Vibration: Millikan-White (1963) + Park high-T collision-limiting correction.
Chemistry: Park (1990) dissociation rates, molecular collision partners.
"""

import math

from .constants import SPECIES, K_B, AMU, ATM, N_A


def _tau_mw(T, p, theta, mu_amu):
    """Millikan-White relaxation time [s]. p in Pa."""
    A = 1.16e-3 * math.sqrt(mu_amu) * theta ** (4.0 / 3.0)
    B = 0.015 * mu_amu ** 0.25
    return math.exp(A * (T ** (-1.0 / 3.0) - B) - 18.42) * ATM / p


def _tau_park(T, n, m_kg):
    """Park high-temperature correction [s]."""
    sigma = 1e-21 * (50000.0 / T) ** 2
    cbar = math.sqrt(8.0 * K_B * T / (math.pi * m_kg))
    return 1.0 / (n * sigma * cbar)


def tau_vib(species, T, p):
    """Mixture-averaged vibrational relaxation time for `species` [s]."""
    s = SPECIES[species]
    m_s = s["M"] / N_A
    n = p / (K_B * T)
    inv = 0.0
    for r in SPECIES.values():
        mu = (s["M"] * r["M"]) / (s["M"] + r["M"]) * 1e3  # amu
        inv += r["X"] / _tau_mw(T, p, s["theta_v"], mu)
    return 1.0 / inv + _tau_park(T, n, m_s)


# Park 1990, partner = molecule: k = C T^n exp(-theta_d/T) [cm^3/(mol s)]
_PARK_DISS = {
    "O2": (2.0e21, -1.5, 59500.0),
    "N2": (7.0e21, -1.6, 113200.0),
}


def tau_diss(species, T, p):
    """Dissociation time 1/(k_f [M]) [s], T as controlling temperature."""
    C, nexp, th = _PARK_DISS[species]
    kf = C * T ** nexp * math.exp(-th / T)          # cm^3/(mol s)
    conc = p / (8.314462618 * T) * 1e-6             # mol/cm^3
    return 1.0 / (kf * conc) if kf > 0 else math.inf


def flow_time(L, u):
    return L / u


def damkohler(t_flow, tau):
    return t_flow / tau if tau > 0 else math.inf
