"""US Standard Atmosphere 1976, 0-86 km (geometric)."""

import math

from .constants import R_MIX

_G0 = 9.80665
_R_EARTH = 6356766.0  # m, US76 effective radius

# (geopotential base [m], lapse rate [K/m])
_LAYERS = [
    (0.0, -6.5e-3),
    (11000.0, 0.0),
    (20000.0, 1.0e-3),
    (32000.0, 2.8e-3),
    (47000.0, 0.0),
    (51000.0, -2.8e-3),
    (71000.0, -2.0e-3),
    (84852.0, None),
]


def _base_states():
    """Precompute T, p at each layer base."""
    T, p = 288.15, 101325.0
    states = [(T, p)]
    for (h0, lam), (h1, _) in zip(_LAYERS[:-1], _LAYERS[1:]):
        dh = h1 - h0
        if lam == 0.0:
            p = p * math.exp(-_G0 * dh / (R_MIX * T))
            T1 = T
        else:
            T1 = T + lam * dh
            p = p * (T1 / T) ** (-_G0 / (R_MIX * lam))
        T = T1
        states.append((T, p))
    return states


_BASES = _base_states()


def us76(alt_m):
    """Return (T [K], p [Pa], rho [kg/m^3]) at geometric altitude alt_m."""
    if not 0.0 <= alt_m <= 86000.0:
        raise ValueError(f"altitude {alt_m/1e3:.1f} km outside US76 range 0-86 km")
    H = _R_EARTH * alt_m / (_R_EARTH + alt_m)  # geopotential
    for i in range(len(_LAYERS) - 1):
        h0, lam = _LAYERS[i]
        if H <= _LAYERS[i + 1][0] or i == len(_LAYERS) - 2:
            Tb, pb = _BASES[i]
            dh = H - h0
            if lam == 0.0:
                T = Tb
                p = pb * math.exp(-_G0 * dh / (R_MIX * Tb))
            else:
                T = Tb + lam * dh
                p = pb * (T / Tb) ** (-_G0 / (R_MIX * lam))
            return T, p, p / (R_MIX * T)
    raise RuntimeError("unreachable")
