"""Hypersonic thermochemical regime estimator.

Estimates total enthalpy, post-shock state, and vibrational / chemical
Damkohler numbers to suggest perfect-gas, 1-T, or 2-T modeling.
"""

from .classify import evaluate, Result  # noqa: F401
