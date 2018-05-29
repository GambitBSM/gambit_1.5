#!/usr/bin/env python

from __future__ import division
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

from read_RHN_chains import *

OUTPATH = '/home/ubuntu/'

def show_survival_fraction(rhn, sigma = 2, exclude = []):
    if len(exclude) > 0:
        lnL = rhn.lnL[rhn.valid]*0
        for key in rhn.lnL_partial:
            if not any([tag in key for tag in exclude]):
                lnL += rhn.lnL_partial[key][rhn.valid]
        lnLmax = lnL.max()
        surv = lnLmax-lnL < 0.5*sigma**2
        print "All - excluded:", surv.sum()/len(lnL)
        print "(excluding:", exclude, ")"
        return

    for name in rhn.lnL_partial:
        lnL = rhn.lnL_partial[name][rhn.valid]
        lnLmax = lnL.max()
        surv = lnLmax-lnL < 0.5*sigma**2
        print name, surv.sum()/len(lnL)

    lnL = rhn.lnL[rhn.valid]
    lnLmax = lnL.max()
    surv = lnLmax-lnL < 0.5*sigma**2
    print "Total (LogLike):", surv.sum()/len(lnL)


def show_lnL_hist(rhn):
    def f(rhn, name):
        print name
        if name == 'total':
            lnL = rhn.lnL[rhn.valid]
        else:
            lnL = rhn.lnL_partial[name][rhn.valid]
        lnL -= lnL.max()
        plt.clf()

        plt.subplot(211)
        plt.xlim([-50, 1])
        frac = (lnL > -100).sum()/len(lnL)
        plt.title(name + " %.2f"%frac)
        plt.hist(lnL, log=True, bins = 200, range = [-50, 1])

        plt.subplot(212)
        plt.xlabel("renormalized lnL")
        frac = (lnL > -100).sum()/len(lnL)
        plt.title(name + " all")
        plt.hist(lnL, log=True, bins = 200)

        plt.savefig(OUTPATH + 'hist_NH_'+name+'.pdf')
    for name in rhn.lnL_partial:
        f(rhn, name)
    f(rhn, 'total')

def show_md21(rhn):
    y = rhn.lnL_partial['lnL_md21'][rhn.valid]
    x = rhn.md21[rhn.valid]*1e18
    plt.clf()
    plt.scatter(x, y, marker = '.')
    plt.xlabel('md21 [eV^2]')
    plt.ylabel('lnL_md21')
    plt.xlim([0, 0.0015])
    plt.savefig(OUTPATH + 'lnL_md21.png')

def show_lnL_inv_Z_width(rhn):
    lnL = rhn.lnL_partial['lnL_Z_inv_width'][rhn.valid]
    U = rhn.U[rhn.valid]
    M1 = rhn.M1[rhn.valid]
    M2 = rhn.M2[rhn.valid]
    M3 = rhn.M3[rhn.valid]

    plt.subplot(221)
    plt.scatter(np.log10(U), lnL)
    plt.xlabel("log10(U)")
    plt.ylabel("lnL_Z_inv_width")

    plt.subplot(222)
    plt.scatter(np.log10(M1), lnL)
    plt.xlabel("log10(M1)")
    plt.ylabel("lnL_Z_inv_width")

    plt.subplot(223)
    plt.scatter(np.log10(M2), lnL)
    plt.xlabel("log10(M2)")
    plt.ylabel("lnL_Z_inv_width")

    plt.subplot(224)
    plt.scatter(np.log10(M3), lnL)
    plt.xlabel("log10(M3)")
    plt.ylabel("lnL_Z_inv_width")

    plt.savefig(OUTPATH+"lnL_Z_inv_width.png")

def check_sum(rhn, exclude = []):
    lnL1 = rhn.lnL[rhn.valid]
    lnL2 = lnL1*0
    for key in rhn.lnL_partial:
        if not any([tag in key for tag in exclude]):
            lnL2 += rhn.lnL_partial[key][rhn.valid]
    print lnL1
    print lnL2

if __name__ == "__main__":
    rhn = RHN_Chain('/home/ubuntu/data/chains/RHN_NH_diff_123.hdf5', 'diff', print_keys = False)
    #check_sum(rhn, exclude = ['inv'])
    #show_lnL_inv_Z_width(rhn)
    #show_md21(rhn)
    show_survival_fraction(rhn)
    show_survival_fraction(rhn, exclude = ['inv', 'LUV'])
    show_survival_fraction(rhn, exclude = ['inv', 'LUV', 'md21', 
        'theta13', 'md3l', 'deltaCP', 'theta12', 'theta23'])
    #show_lnL_hist(rhn)
