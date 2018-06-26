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

def show_high_couplings(rhn):
    mask = rhn.lnL > rhn.lnL.max()-2
    M = [rhn.M1, rhn.M2, rhn.M3]
    U = [
        [rhn.U1, rhn.Ue1, rhn.Um1, rhn.Ut1],
        [rhn.U2, rhn.Ue2, rhn.Um2, rhn.Ut2],
        [rhn.U3, rhn.Ue3, rhn.Um3, rhn.Ut3],
        ]
    quit()

def show_lnL_relevance(rhn, i = 1, I = 1, tag = 'TAG'):
    M = [rhn.M1, rhn.M2, rhn.M3]
    U = [
        [rhn.U1, rhn.Ue1, rhn.Um1, rhn.Ut1],
        [rhn.U2, rhn.Ue2, rhn.Um2, rhn.Ut2],
        [rhn.U3, rhn.Ue3, rhn.Um3, rhn.Ut3],
        ]

    mask = rhn.lnL > rhn.lnL.max()-2

    def f(rhn, name):
        print name
        plt.clf()
        lnL = rhn.lnL[mask]

        if name == 'total':
            w = lnL
            indices = np.argsort(w)
        else:
            lnL_partial = rhn.lnL_partial[name][mask]
            imax = lnL.argmax()
            x = lnL[imax] - lnL
            y = lnL_partial[imax] - lnL_partial
            w = y/(x+0.1)
            indices = np.argsort(-w)
        plt.scatter(np.log10(M[I-1])[mask][indices],
                np.log10(U[I-1][i])[mask][indices], c = w[indices],
                rasterized=True, marker='.', edgecolors = None, linewidths = 0)
        plt.colorbar()
        plt.xlabel("log10(M%i/GeV)"%I)
        plt.ylabel("log10(U%i%i)"%(i,I))

        print 'save...'
        plt.savefig(OUTPATH + 'hist_NH_scatter_U%i%i'%(i,I)+"_"+name+'_'+tag+'.pdf')
        print '...done'

    for name in rhn.lnL_partial:
        f(rhn, name)
    #f(rhn, 'total')

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
    print "mNu..."
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
    lnL = rhn.lnL

    def plot(mask, filename):
        #print (mNu2**2-mNu1**2).min()
        #print (mNu2**2-mNu1**2).max()
        #print (mNu3**2-mNu1**2).min()
        #print (mNu3**2-mNu1**2).max()
        #print (mNu3**2-mNu2**2).min()
        #print (mNu3**2-mNu2**2).max()

        plt.clf()
        plt.subplot(221)
        plt.scatter(np.log10(mNu_1st)[mask], np.log10(mNu_2nd)[mask], marker='.', alpha = 0.01, rasterized = True)
        plt.scatter(np.log10(mNu_1st)[mask], np.log10(mNu_3rd)[mask], marker='.', alpha = 0.01, rasterized = True)

        plt.subplot(222)
        plt.hist(np.log10(mNu1)[mask], bins = 200, log=True)
        plt.hist(np.log10(mNu2)[mask], bins = 200, log=True)
        plt.hist(np.log10(mNu3)[mask], bins = 200, log=True)

        plt.subplot(223)
        plt.scatter(np.log10(mNu_1st)[mask], np.log10(md21)[mask], marker='.', alpha = 0.01, rasterized = True)
        plt.scatter(np.log10(mNu_1st)[mask], np.log10(abs(md3l))[mask], marker='.', alpha = 0.01, rasterized = True)

        plt.subplot(224)
        plt.hist(np.log10(md21)[mask], bins = 200, log=True)
        plt.hist(np.log10(abs(md3l))[mask], bins = 200, log=True)

        plt.savefig(OUTPATH+filename, dpi = 200)

    mask = (lnL.max() - lnL > 2) & (md31 > 0)
    plot(mask, 'mNu_NH.pdf')
    mask = (lnL.max() - lnL > 2) & (md31 < 0)
    plot(mask, 'mNu_IH.pdf')

def show_Rorder(rhn):
    print "Rorder..."
    plt.clf()
    plt.hist(rhn.Rorder, range = [0, 6], bins = 6)
    plt.savefig(OUTPATH+"Rorder.pdf", dpi = 200)

def show_U_vs_M(rhn, tag = "TAG"):
    print "U_vs_M..."
    lnL = rhn.lnL
    mask = lnL.max()-lnL < 2
    mask_ft = get_protected(rhn, epsilon = 1e-3, eta = 1e-3, mbbK = 1e-3)
    mask2 = mask & mask_ft

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
        plt.scatter(np.log10(M[I-1])[mask2], np.log10(U[I-1][0])[mask2], marker = '.',
                rasterized = True, color='g')
        plt.xlim([-1, 3.0])
        plt.ylim([-10, -1])
        plt.ylabel("U%i"%I)
        plt.subplot(222)
        plt.scatter(np.log10(M[I-1])[mask], np.log10(U[I-1][1])[mask], marker = '.',
                rasterized = True)
        plt.scatter(np.log10(M[I-1])[mask2], np.log10(U[I-1][1])[mask2], marker = '.',
                rasterized = True, color='g')
        plt.xlim([-1, 3.0])
        plt.ylim([-10, -1])
        plt.ylabel("Ue%i"%I)
        plt.subplot(223)
        plt.scatter(np.log10(M[I-1])[mask], np.log10(U[I-1][2])[mask], marker = '.',
                rasterized = True)
        plt.scatter(np.log10(M[I-1])[mask2], np.log10(U[I-1][2])[mask2], marker = '.',
                rasterized = True, color='g')
        plt.xlim([-1, 3.0])
        plt.ylim([-10, -1])
        plt.ylabel("Um%i"%I)
        plt.subplot(224)
        plt.scatter(np.log10(M[I-1])[mask], np.log10(U[I-1][3])[mask], marker = '.',
                rasterized = True)
        plt.scatter(np.log10(M[I-1])[mask2], np.log10(U[I-1][3])[mask2], marker = '.',
                rasterized = True, color='g')
        plt.ylabel("Ut%i"%I)
        plt.xlim([-1, 3.0])
        plt.ylim([-10, -1])
        plt.tight_layout(pad=0.3)
        plt.savefig(OUTPATH+"U_vs_M%i_%s.pdf"%(I,tag))

def show_lnL_mbb(rhn):
    lnL0 = -rhn.lnL_partial['lnL_0nubb']*2
    lnL1 = -rhn.lnL_partial['lnL_mbb_0nubb_KamLAND_Zen']*2
    lnL2 = -rhn.lnL_partial['lnL_mbb_0nubb_GERDA']*2
    lnL0 -= lnL0.min()
    lnL1 -= lnL1.min()
    lnL2 -= lnL2.min()
    plt.scatter(np.log10(lnL0), np.log10(lnL1+lnL2), marker='.', color='r',
            label="KamLAND_Zen", rasterized = True)
    plt.scatter(np.log10(lnL0), np.log10(lnL2), marker='.', color='g',
            label="GERDA", rasterized = True)
    plt.xlim([-3, 3])
    plt.ylim([-3, 3])
    plt.xlabel("log10(chi2_old)")
    plt.ylabel("log10(chi2_new)")
    plt.axhline(0)
    plt.axvline(0)
    plt.legend()
    plt.savefig(OUTPATH+"lnL_0nubb.pdf")

def show_0nubb_impact(rhn):
    lnL0 = -rhn.lnL_partial['lnL_0nubb']*2
    lnL1 = -rhn.lnL_partial['lnL_mbb_0nubb_KamLAND_Zen']*2
    lnL2 = -rhn.lnL_partial['lnL_mbb_0nubb_GERDA']*2
    lnL0 -= lnL0.min()
    lnL1 -= lnL1.min()
    lnL2 -= lnL2.min()

    plt.scatter(np.log10(rhn.U), np.log10(lnL0), rasterized = True, marker='.',
            color='r')
    plt.scatter(np.log10(rhn.U), np.log10(lnL1), rasterized = True, marker='.',
            color='g')
    plt.axhline(0)
    plt.axvline(-9)
    plt.axvline(-8)

    plt.savefig(OUTPATH+"U_vs_M_0nubb.pdf")

def show_phases(rhn):
    plt.scatter(np.angle(rhn.ue1), np.angle(rhn.um1), rasterized = True, marker='.')
    plt.savefig(OUTPATH+"angles.pdf")

def show_mbb(rhn):
    print "mbb..."
    plt.clf()
    lnL = rhn.lnL
    U = rhn.U
    md31 = rhn.md31
    mask = (lnL.max() - lnL < 2)

    # Exclude non-protected points
    mask_ft = get_protected(rhn, epsilon = 1e-2, eta = 1e-0, mbbK = 1e-4)
    mask &= mask_ft

    mask1 = mask & (md31 > 0)
    print "Number of points NH:", mask1.sum()
    mask2 = mask & (md31 < 0)
    print "Number of points IH:", mask2.sum()
    mbb = rhn.mbb*1e9
    mMin = rhn.mMin*1e9
    plt.scatter(np.log10(mMin), np.log10(mbb), rasterized = True, marker='.', color='0.5')
    plt.scatter(np.log10(mMin)[mask1], np.log10(mbb)[mask1], rasterized = True, marker='.', color='g')
    plt.scatter(np.log10(mMin)[mask2], np.log10(mbb)[mask2], rasterized = True, marker='.', color='r')
    plt.ylim([-4, 0])
    plt.xlim([-6, -1])
    plt.xlabel("log10(m_light [eV])")
    plt.ylabel("log10(mbb [eV])")
    plt.savefig(OUTPATH+"mbb.pdf")

def get_protected(rhn, IJK = None, epsilon = np.inf, eta = np.inf, mbbK = np.inf):
    M = [rhn.M1, rhn.M2, rhn.M3]
    ue = [rhn.ue1, rhn.ue2, rhn.ue3]

    def f(I, J, K):
        dM = abs(M[I] - M[J])
        m1 = dM < epsilon*(M[I]+M[J])/2
        m1 &= M[I] > 10.
        deta = abs(ue[I]**2 + ue[J]**2)/(abs(ue[I])**2+abs(ue[J])**2)
        m2 = deta < eta
        p2 = 0.1**2  # 100 MeV  (ad hoc)
        mbb_I = p2*abs(ue[I])**2/M[I]*1e9  # eV
        mbb_J = p2*abs(ue[J])**2/M[J]*1e9  # eV
        mbb_K = p2*abs(ue[K])**2/M[K]*1e9  # eV
    #    mbb += p2*abs(ue[I])**2/M[I]*1e9  # eV
    #    mbb += p2*abs(ue[J])**2/M[J]*1e9  # eV
        #m3 = (mbb_K < mbbK) & (mbb_I > 0.01)
        m3 = mbb_K < mbbK  # eV
    #    print mbb[m1 & m2].min()
    #    print mbb[m1 & m2].max()
    #    print len(mbb[m1 & m2])
    #    print (mbb[m1 & m2] > mbbK).sum()
    #    quit()
        print "First cut:", m1.sum()/len(m1)
        print "Second cut:", (m1&m2).sum()/m1.sum()
        print "Third cut:", (m1&m2&m3).sum()/(m1&m2).sum()
        mtot = m1 & m2 & m3
        return mtot

    if IJK is None:
        m123 = f(0, 1, 2)
        m231 = f(1, 2, 0)
        m312 = f(2, 0, 1)
        mtot = m123 | m231 | m312
    else:
        if IJK not in [123, 231, 312, 321, 213, 132]:
            raise KeyError("IJK invalid.")
        IJK = str(IJK)
        I, J, K = int(IJK[0])-1, int(IJK[1])-1, int(IJK[2])-1
        mtot = f(I, J, K)

    return mtot

def print_finetuning_counts(rhn):
    print get_protected(rhn, IJK = 123, epsilon = 1e-2).sum()
    print get_protected(rhn, IJK = 231, epsilon = 1e-2).sum()
    print get_protected(rhn, IJK = 312, epsilon = 1e-2).sum()
    print get_protected(rhn, epsilon = 1e-2).sum()
    print get_protected(rhn, IJK = 123, epsilon = 1e-2, eta = 1e-2).sum()
    print get_protected(rhn, IJK = 231, epsilon = 1e-2, eta = 1e-2).sum()
    print get_protected(rhn, IJK = 312, epsilon = 1e-2, eta = 1e-2).sum()
    print get_protected(rhn, epsilon = 1e-2, eta = 1e-2).sum()
    print get_protected(rhn, IJK = 123, epsilon = 1e-2, eta = 1e-2, mbbK = 0.1).sum()
    print get_protected(rhn, IJK = 231, epsilon = 1e-2, eta = 1e-2, mbbK = 0.1).sum()
    print get_protected(rhn, IJK = 312, epsilon = 1e-2, eta = 1e-2, mbbK = 0.1).sum()
    print get_protected(rhn, epsilon = 1e-2, eta = 1e-2, mbbK = 0.1).sum()

def triangle(rhn):
    print "triangle..."
    x = rhn.Ue1/rhn.U1
    y = rhn.Um1/rhn.U1
    mMin = rhn.mMin * 1e9
    lnL = rhn.lnL
    md31 = rhn.md31

    plt.subplot(121)
    mask = (lnL.max() - lnL < 2) & (mMin < 0.01) & (md31 > 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='k', marker='.')
    mask = (lnL.max() - lnL < 2) & (mMin < 0.001) & (md31 > 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='r', marker='.')
    mask = (lnL.max() - lnL < 2) & (mMin < 0.0001) & (md31 > 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='g', marker='.')
    mask = (lnL.max() - lnL < 2) & (mMin < 0.00001) & (md31 > 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='y', marker='.')
    plt.xlabel("Ue1/U1")
    plt.ylabel("Um1/U1")

    plt.subplot(122)
    mask = (lnL.max() - lnL < 2) & (mMin < 0.01) & (md31 < 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='k', marker='.')
    mask = (lnL.max() - lnL < 2) & (mMin < 0.001) & (md31 < 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='r', marker='.')
    mask = (lnL.max() - lnL < 2) & (mMin < 0.0001) & (md31 < 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='g', marker='.')
    mask = (lnL.max() - lnL < 2) & (mMin < 0.00001) & (md31 < 0)
    plt.scatter(x[mask], y[mask], rasterized = True, color='y', marker='.')
    plt.xlabel("Ue1/U1")
    plt.ylabel("Um1/U1")

    plt.savefig(OUTPATH+"triangle.pdf")

if __name__ == "__main__":
    rhn = RHN_Chain('/home/ubuntu/data/RHN_diff_NH_cs18.hdf5', MODEL = 'diff',
            print_keys = False, renormalize = False)
    #triangle(rhn)
    #show_mbb(rhn)
    #show_Rorder(rhn)
    show_U_vs_M(rhn, tag = 'cs18')
    #show_survival_fraction(rhn)
    #show_neutrino_masses(rhn)
    #show_high_couplings(rhn)

    #print_finetuning_counts(rhn)
    #show_phases(rhn)
    #show_lnL_mbb(rhn)
    #show_0nubb_impact(rhn)
    #check_sum(rhn, exclude = ['inv'])
    #show_lnL_inv_Z_width(rhn)
    #show_md21(rhn)
    #show_survival_fraction(rhn, exclude = ['inv', 'LUV'])
    #show_survival_fraction(rhn, exclude = ['inv', 'LUV', 'md21', 
    #    'theta13', 'md3l', 'deltaCP', 'theta12', 'theta23'])
    #show_lnL_hist(rhn)c
    #show_lnL_relevance(rhn, i = 3, I = 1, tag = 'cs15')
