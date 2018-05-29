"""
Profile a chi-squared in a data file with a fractional theory error
"""

import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar


def profile(file_name, frac_error=0.1, min_=0., max_=1.):
    """
    Profile a chi-squared in a data file with a fractional theory error
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

    FILE_NAME = sys.argv[1]

    try:
        FRAC_ERROR = float(sys.argv[2])
    except IndexError:
        FRAC_ERROR = 0.1

    try:
        MIN = float(sys.argv[3])
    except IndexError:
        MIN = 0.

    try:
        MAX = float(sys.argv[4])
    except IndexError:
        MAX = 1.

    for p, c in profile(FILE_NAME, FRAC_ERROR, MIN, MAX):
        print p, c
