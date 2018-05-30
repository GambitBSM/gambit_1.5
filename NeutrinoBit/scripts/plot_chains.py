#!/usr/bin/env python

from __future__ import division
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

from read_RHN_chains import *

OUTPATH = '/home/cweniger/'

def show_survival_fraction(rhn, sigma = 2, exclude = []):
    if len(exclude) > 0:
        lnL = rhn.lnL*0
        for key in rhn.lnL_partial:
            if not any([tag in key for tag in exclude]):
                lnL += rhn.lnL_partial[key]
        lnLmax = lnL.max()
        surv = lnLmax-lnL < 0.5*sigma**2
        print "All - excluded:", surv.sum()/len(lnL)
        print "(excluding:", exclude, ")"
        return

    for name in rhn.lnL_partial:
        lnL = rhn.lnL_partial[name]
        lnLmax = lnL.max()
        surv = lnLmax-lnL < 0.5*sigma**2
        print name, surv.sum()/len(lnL)

    lnL = rhn.lnL
    lnLmax = lnL.max()
    surv = lnLmax-lnL < 0.5*sigma**2
    print "Total (LogLike):", surv.sum()/len(lnL)


def show_lnL_hist(rhn):
    def f(rhn, name):
        print name
        if name == 'total':
            lnL = rhn.lnL
        else:
            lnL = rhn.lnL_partial[name]
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
    y = rhn.lnL_partial['lnL_md21']
    x = rhn.md21*1e18
    plt.clf()
    plt.scatter(x, y, marker = '.')
    plt.xlabel('md21 [eV^2]')
    plt.ylabel('lnL_md21')
    plt.xlim([0, 0.0015])
    plt.savefig(OUTPATH + 'lnL_md21.png')

def show_lnL_inv_Z_width(rhn):
    lnL = rhn.lnL_partial['lnL_Z_inv_width']
    U = rhn.U
    M1 = rhn.M1
    M2 = rhn.M2
    M3 = rhn.M3

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
    lnL1 = rhn.lnL
    lnL2 = lnL1*0
    for key in rhn.lnL_partial:
        if not any([tag in key for tag in exclude]):
            lnL2 += rhn.lnL_partial[key]
    print lnL1
    print lnL2

def show_neutrino_masses(rhn):
    mNu1 = rhn.mNu1*1e9
    mNu2 = rhn.mNu2*1e9
    mNu3 = rhn.mNu3*1e9
    md21 = rhn.md21*1e18
    md32 = rhn.md32*1e18
    md31 = rhn.md31*1e18
    md3l = np.where(md32 > 0, md32, md31)
    mNu_1st = np.minimum(mNu1, mNu3)
    mNu_2nd = np.where(mNu1 < mNu3, mNu2, mNu1)
    mNu_3rd = np.where(mNu1 < mNu3, mNu3, mNu2)

    #print (mNu2**2-mNu1**2).min()
    #print (mNu2**2-mNu1**2).max()
    #print (mNu3**2-mNu1**2).min()
    #print (mNu3**2-mNu1**2).max()
    #print (mNu3**2-mNu2**2).min()
    #print (mNu3**2-mNu2**2).max()

    plt.clf()
    plt.subplot(221)
    plt.scatter(np.log10(mNu_1st), np.log10(mNu_2nd), marker='.', alpha = 0.01, rasterized = True)
    plt.scatter(np.log10(mNu_1st), np.log10(mNu_3rd), marker='.', alpha = 0.01, rasterized = True)

    plt.subplot(222)
    plt.hist(np.log10(mNu1), bins = 200, log=True)
    plt.hist(np.log10(mNu2), bins = 200, log=True)
    plt.hist(np.log10(mNu3), bins = 200, log=True)

    plt.subplot(223)
    plt.scatter(np.log10(mNu_1st), np.log10(md21), marker='.', alpha = 0.01, rasterized = True)
    plt.scatter(np.log10(mNu_1st), np.log10(abs(md3l)), marker='.', alpha = 0.01, rasterized = True)

    plt.subplot(224)
    plt.hist(np.log10(md21), bins = 200, log=True)
    plt.hist(np.log10(abs(md3l)), bins = 200, log=True)

    plt.savefig(OUTPATH+"mNu.pdf", dpi = 200)

def show_Rorder(rhn):
    plt.clf()
    plt.subplot(211)
    plt.scatter(rhn.Rorder, np.log10(rhn.Ue1), marker='.', rasterized = True)
    plt.subplot(212)
    plt.hist(rhn.Rorder, range = [0, 6], bins = 6)
    plt.savefig(OUTPATH+"Rorder.pdf", dpi = 200)

def show_U_vs_M(rhn):
    lnL = rhn.lnL
    mask = lnL.max()-lnL < 2

    M = [rhn.M1, rhn.M2, rhn.M3]
    U = [
            [rhn.U1, rhn.Ue1, rhn.Um1, rhn.Ut1],
            [rhn.U2, rhn.Ue2, rhn.Um2, rhn.Ut2],
            [rhn.U3, rhn.Ue3, rhn.Um3, rhn.Ut3],
            ]

    for I in [1, 2, 3]:
        print "Generating U_vs_M%i"%I
        plt.clf()
        plt.subplot(221)
        plt.scatter(np.log10(M[I-1])[mask], np.log10(U[I-1][0])[mask], marker = '.',
                rasterized = True)
        plt.xlim([-1, 3.0])
        plt.ylim([-17, -3])
        plt.ylabel("U%i"%I)
        plt.subplot(222)
        plt.scatter(np.log10(M[I-1])[mask], np.log10(U[I-1][1])[mask], marker = '.',
                rasterized = True)
        plt.xlim([-1, 3.0])
        plt.ylim([-17, -3])
        plt.ylabel("Ue%i"%I)
        plt.subplot(223)
        plt.scatter(np.log10(M[I-1])[mask], np.log10(U[I-1][2])[mask], marker = '.',
                rasterized = True)
        plt.xlim([-1, 3.0])
        plt.ylim([-17, -3])
        plt.ylabel("Um%i"%I)
        plt.subplot(224)
        plt.scatter(np.log10(M[I-1])[mask], np.log10(U[I-1][3])[mask], marker = '.',
                rasterized = True)
        plt.ylabel("Ut%i"%I)
        plt.xlim([-1, 3.0])
        plt.ylim([-17, -3])
        plt.tight_layout(pad=0.3)
        plt.savefig(OUTPATH+"U_vs_M%i.pdf"%I)

def show_lnL_mbb(rhn):
    lnL0 = rhn.lnL_partial['lnL_0nubb']
    lnL1 = rhn.lnL_partial['lnL_mbb_0nubb_KamLAND_Zen']
    lnL2 = rhn.lnL_partial['lnL_mbb_0nubb_GERDA']
    lnL0 -= lnL0.max()
    lnL1 -= lnL1.max()
    lnL2 -= lnL2.max()

    plt.scatter(-lnL0, -lnL1, marker='.', label="KamLAND_Zen", rasterized = True)
    plt.scatter(-lnL0, -lnL2, marker='.', label="GERDA", rasterized = True)
    plt.xlim([0, 100])
    plt.ylim([0, 100])
    plt.legend()
    plt.savefig(OUTPATH+"lnL_0nubb.pdf")

if __name__ == "__main__":
    #rhn = RHN_Chain('/home/cweniger/hdf5_29_05_2018/RHN_diff_NH_123_1e-5.hdf5', print_keys = False)
    rhn = RHN_Chain('/home/cweniger/hdf5_29_05_2018/RHN_diff_NH_1e-9_old0nubb.hdf5', print_keys = False)
    #show_Rorder(rhn)
    #show_neutrino_masses(rhn)
    show_lnL_mbb(rhn)
    #check_sum(rhn, exclude = ['inv'])
    #show_lnL_inv_Z_width(rhn)
    #show_md21(rhn)
    #show_U_vs_M(rhn)
    #show_survival_fraction(rhn)
    #show_survival_fraction(rhn, exclude = ['inv', 'LUV'])
    #show_survival_fraction(rhn, exclude = ['inv', 'LUV', 'md21', 
    #    'theta13', 'md3l', 'deltaCP', 'theta12', 'theta23'])
    #show_lnL_hist(rhn)
