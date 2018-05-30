#!/usr/bin/env python

from __future__ import division
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt


def symmetrize(infile, outfile):
    print "WARNING: Only symmetrizes RHN masses and couplings and phases."

    root = h5py.File(infile)
    group = root["RHN"]
    keys = group.keys()

    def get_data(tag):
        return np.array(group[tag])

    vector_tags = [
            lambda I: '#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::M_%i'%I,
            lambda I: '#Ue%i @NeutrinoBit::Ue%i'%(I,I),
            lambda I: '#Ue%i_phase @NeutrinoBit::Ue%i_phase'%(I,I),
            lambda I: '#Um%i @NeutrinoBit::Um%i'%(I,I),
            lambda I: '#Um%i_phase @NeutrinoBit::Um%i_phase'%(I,I),
            lambda I: '#Ut%i @NeutrinoBit::Ut%i'%(I,I),
            lambda I: '#Ut%i_phase @NeutrinoBit::Ut%i_phase'%(I,I),
            ]

    key_to_vector_tag = dict()
    key_to_vector_index = dict()
    for i, key in enumerate(keys):
        key_to_vector_tag[key] = None
        key_to_vector_index[key] = 0
        for vector_tag in vector_tags:
            for I in [1, 2, 3]:
                if key == vector_tag(I):
                    key_to_vector_tag[key] = vector_tag
                    key_to_vector_index[key] = I
    Iperm = [
            [1, 2, 3, 1, 2, 3],
            [2, 3, 1, 3, 1, 2],
            [3, 1, 2, 2, 3, 1],
            ]

    hf = h5py.File(outfile, 'w')
    outgroup = hf.create_group("RHN")
    for key in keys:
        keyperm = [key, key, key, key, key, key]
        if key_to_vector_tag[key] is not None:
            I = key_to_vector_index[key]
            keyperm = [key_to_vector_tag[key](J) for J in Iperm[I-1]]
            print key, "......................................(vector)"
        else:
            print key
        data = np.concatenate([np.array(group[keyp]) for keyp in keyperm])
        outgroup.create_dataset(key, data = data)
    hf.close()
    root.close()

if __name__ == "__main__":
    symmetrize('/home/cweniger/hdf5_29_05_2018/RHN_diff_IH_1e-8.hdf5', 'test.hdf5')
