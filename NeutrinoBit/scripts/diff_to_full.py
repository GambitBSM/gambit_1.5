#!/usr/bin/env python

from __future__ import division
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

def diff_to_full(infile, outfile):
    root = h5py.File(infile)
    group = root["RHN"]
    keys = group.keys()

    M1 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_1'])
    #M3 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_3'])
    dM21 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::delta_M_2'])

    M2 = M1 + dM21
    M3 = M1 + dM31
    #M3 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_3'])

    key_M2 = '#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_2'
    key_M3 = '#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_2'

    hf = h5py.File(outfile, 'w')
    outgroup = hf.create_group("RHN")
    for key in keys:
        data = group[key]
        outgroup.create_dataset(key, data = data)
    outgroup.create_dataset(key_M2, data = M2)
    hf.close()
    root.close()

if __name__ == "__main__":
    diff_to_full('/home/cweniger/hdf5_29_05_2018/RHN_diff_IH_1e-8.hdf5', 'test.hdf5')
