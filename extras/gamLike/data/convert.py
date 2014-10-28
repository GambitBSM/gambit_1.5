#!/usr/bin/env python

from __future__ import division
from numpy import *
import pylab as plt
from scipy import linalg

data = loadtxt('covariance.dat')
out = zeros((24, 27))
out[:,0] = data[:,0]
out[:,1] = data[:,1]
out[:,2] = data[:,2]
sigma = data[:,-24:]
out[:,3:] = linalg.inv(sigma)

emeans = (out[:,0] * out[:,1])**0.5
de = (out[:,1] - out[:,0])

out[:,2] *= emeans**-2 * de
a, b = meshgrid(emeans**-2 * de, emeans**-2 * de)
out[:,3:] /= a*b

savetxt('like_GC_Calore.txt', out)
