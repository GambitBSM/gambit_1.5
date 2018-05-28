"""
Convolve chi-squared in a data file with a fractional theory error
"""

import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
from scipy.integrate import quad


def convolve(file_name, frac_error=0.1, min_=0., max_=1.):
    """
    Convolve chi-squared in a data file with a fractional theory error
    """
    # Unpack data
    param, chi_squared = np.loadtxt(file_name, unpack=True)

    # Interpolate likelihood function
    like = interp1d(param, np.exp(-0.5 * chi_squared),
        kind='linear', bounds_error=False, fill_value="extrapolate")

    # Make prior for true prediction
    gaussian = lambda x, mu: norm.pdf(x, mu, frac_error * mu)

    # Make functions for convolution
    integrand = lambda x, mu: like(x) * gaussian(x, mu)
    convolved = lambda mu: -2. * np.log(
        quad(integrand, min_, max_, args=(mu))[0]
        / quad(gaussian, min_, max_, args=(mu))[0])

    # Perform convolution
    return [(p, convolved(p)) for p in param]

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

    for p, c in convolve(FILE_NAME, FRAC_ERROR, MIN, MAX):
        print p, c
