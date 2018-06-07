#!/usr/bin/env python
"""
Profile a chi-squared in a data file with a fractional theory error
===================================================================

python profile_theory.py <file> <frac-error> <min> <max>

prints (parameter, profiled chi-squared) from a file containing
(parameter, chi-squared).
"""

import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar


def profile(file_name, frac_error=0.1, min_=0., max_=1.):
    """
    Profile a chi-squared in a data file with a fractional theory error

    Args:
        file_name (str): Data file with columns (parameter, chi-squared).
        frac_error (float, optional): Fractional theory error on the parameter.
        min_ (float, optional): Minimum value of parameter.
        max_ (float, optional): Maximum value of parameter.

    Returns:
        list(tuples): List of (parameter, profiled chi-squared)
    """
    # Unpack data
    param, chi_squared_ = np.loadtxt(file_name, unpack=True)

    # Interpolate chi-squared function
    chi_squared = interp1d(param, chi_squared_,
        kind='linear', bounds_error=False, fill_value="extrapolate")

    # Make penalty for true prediction
    penalty = lambda x, mu: (x - mu)**2 / (frac_error * mu)**2

    # Make functions for profile
    objective = lambda x, mu: chi_squared(x) + penalty(x, mu)
    profiled = lambda mu: minimize_scalar(objective, bracket=[min_, mu, max_],
      tol=1E-10, args=(mu)).fun

    # Profile
    return [(p, profiled(p)) for p in param]

if __name__ == "__main__":

    try:
        FILE_NAME = sys.argv[1]
    except IndexError:
        raise RuntimeError(__doc__)

    ARGS = map(float, sys.argv[2:])

    for p, c in profile(FILE_NAME, *ARGS):
        print p, c
