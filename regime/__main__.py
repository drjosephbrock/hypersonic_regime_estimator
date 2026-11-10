"""CLI: python -m regime --alt 40 --vel 3000 [--rnose 0.05]

Ranges as start:stop:step, e.g. --alt 20:70:10 --vel 2000:7000:1000
"""

import argparse

import numpy as np

from .classify import evaluate
from . import report


def _rng(s):
    if ":" in s:
        a, b, c = (float(x) for x in s.split(":"))
        return list(np.arange(a, b + 0.5 * c, c))
    return [float(x) for x in s.split(",")]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="regime", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alt", default="40", help="altitude km (value, list, or a:b:step)")
    ap.add_argument("--vel", default="3000", help="velocity m/s (value, list, or a:b:step)")
    ap.add_argument("--mach", help="use Mach instead of --vel (value, list, or a:b:step)")
    ap.add_argument("--rnose", type=float, default=0.1, help="nose radius m (default 0.1)")
    ap.add_argument("--length", type=float, help="override flow length m (default: standoff)")
    ap.add_argument("--T1", type=float, help="freestream T K (e.g. tunnel); with --p1")
    ap.add_argument("--p1", type=float, help="freestream p Pa; with --T1")
    ap.add_argument("--map", metavar="PNG", help="write V-alt regime map to PNG")
    a = ap.parse_args(argv)

    if a.map:
        from .plot import regime_map
        print("wrote", regime_map(np.arange(10, 86, 2), np.arange(500, 8001, 125),
                                  R_nose=a.rnose, fname=a.map))
        return

    alts = _rng(a.alt)
    results = []
    for alt in alts:
        if a.mach:
            from .atmosphere import us76
            from .thermo import gamma
            from .constants import R_MIX
            T1 = a.T1 if a.T1 else us76(alt * 1e3)[0]
            c = (gamma(T1) * R_MIX * T1) ** 0.5
            vels = [m * c for m in _rng(a.mach)]
        else:
            vels = _rng(a.vel)
        for v in vels:
            results.append(evaluate(alt, v, R_nose=a.rnose, L=a.length, T1=a.T1, p1=a.p1))

    print(report.single(results[0]) if len(results) == 1 else report.table(results))


if __name__ == "__main__":
    main()
