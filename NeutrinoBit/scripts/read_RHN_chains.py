#!/usr/bin/env python

from __future__ import division
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

OUTPATH = '/home/cweniger/'

class RHN_Chain(object):
    def __init__(self, INFILE, MODEL):
        print "Reading %s..."%INFILE
        root = h5py.File(INFILE)
        group = root["RHN"]

        #for key in group.keys():
        #    print key
        #quit()

        if MODEL == 'diff':
            self.M1 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_1'])
            self.dM2 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::delta_M_2'])
            self.dM3 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::delta_M_3'])
            self.ReOmega12 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ReOm12'])
            self.ReOmega13 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ReOm13'])
            self.ReOmega23 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ReOm23'])
            self.ImOmega12 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ImOm12'])
            self.ImOmega13 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ImOm13'])
            self.ImOmega23 = np.array(group['#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ImOm23'])
            self.M2 = self.M1 + self.dM2
            self.M3 = self.M1 + self.dM3
        elif MODEL == 'full':
            self.M1 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::M_1'])
            self.M2 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::M_2'])
            self.M3 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::M_3'])
            self.ReOmega12 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ReOm12'])
            self.ReOmega13 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ReOm13'])
            self.ReOmega23 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ReOm23'])
            self.ImOmega12 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ImOm12'])
            self.ImOmega13 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ImOm13'])
            self.ImOmega23 = np.array(group['#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ImOm23'])
            self.dM2 = self.M2 - self.M1
            self.dM3 = self.M3 - self.M1
        else:
            raise KeyError("Model unknown")

        self.Ue1 = np.array(group['#Ue1 @NeutrinoBit::Ue1'])
        self.Ue2 = np.array(group['#Ue2 @NeutrinoBit::Ue2'])
        self.Ue3 = np.array(group['#Ue3 @NeutrinoBit::Ue3'])
        self.Um1 = np.array(group['#Um1 @NeutrinoBit::Um1'])
        self.Um2 = np.array(group['#Um2 @NeutrinoBit::Um2'])
        self.Um3 = np.array(group['#Um3 @NeutrinoBit::Um3'])
        self.Ut1 = np.array(group['#Ut1 @NeutrinoBit::Ut1'])
        self.Ut2 = np.array(group['#Ut2 @NeutrinoBit::Ut2'])
        self.Ut3 = np.array(group['#Ut3 @NeutrinoBit::Ut3'])

        lnL_deltaCP = np.array(group['#deltaCP_lnL @NeutrinoBit::deltaCP_lnL'])
        lnL_l2lgamma = np.array(group['#l2lgamma_lnL @FlavBit::l2lgamma_likelihood'])
        lnL_l2lll = np.array(group['#l2lll_lnL @FlavBit::l2lll_likelihood'])
        lnL_0nubb = np.array(group['#lnL_0nubb @NeutrinoBit::lnL_0nubb'])
        lnL_W_decays = np.array(group['#lnL_W_decays @PrecisionBit::lnL_W_decays_chi2'])
        lnL_W_mass = np.array(group['#lnL_W_mass @PrecisionBit::lnL_W_mass_chi2'])
        lnL_Z_inv_width = np.array(group['#lnL_Z_inv_width @PrecisionBit::lnL_Z_inv_width_chi2'])
        lnL_bbn = np.array(group['#lnL_bbn @NeutrinoBit::lnL_bbn'])
        lnL_atlas_e = np.array(group['#lnLatlase @NeutrinoBit::lnL_atlas_e'])
        lnL_atlas_mu = np.array(group['#lnLatlasmu @NeutrinoBit::lnL_atlas_mu'])
        lnL_charm_e = np.array(group['#lnLcharme @NeutrinoBit::lnL_charm_e'])
        lnL_charm_mu = np.array(group['#lnLcharmmu @NeutrinoBit::lnL_charm_mu'])
        lnL_charm_tau = np.array(group['#lnLcharmtau @NeutrinoBit::lnL_charm_tau'])
        lnL_ckm = np.array(group['#lnLckm_Vusmin @NeutrinoBit::lnL_ckm_Vusmin'])
        lnL_delphi = np.array(group['#lnLdelphi @NeutrinoBit::lnL_delphi'])
        lnL_e949 = np.array(group['#lnLe949 @NeutrinoBit::lnL_e949'])
        lnL_lepuniv = np.array(group['#lnLlepuniv @NeutrinoBit::lnL_lepuniv'])
        lnL_nutev = np.array(group['#lnLnutev @NeutrinoBit::lnL_nutev'])
        lnL_pienu = np.array(group['#lnLpienu @NeutrinoBit::lnL_pienu'])
        lnL_ps191_e = np.array(group['#lnLps191e @NeutrinoBit::lnL_ps191_e'])
        lnL_ps191_mu = np.array(group['#lnLps191mu @NeutrinoBit::lnL_ps191_mu'])
        lnL_md21 = np.array(group['#md21_lnL @NeutrinoBit::md21_lnL'])
        lnL_md3l = np.array(group['#md3l_lnL @NeutrinoBit::md3l_lnL'])
        lnL_mu2e = np.array(group['#mu2e_lnL @FlavBit::mu2e_likelihood'])
        lnL_pert = np.array(group['#perturbativity_lnL @NeutrinoBit::perturbativity_likelihood'])
        lnL_sum_mnu = np.array(group['#sum_mnu_lnL @NeutrinoBit::sum_mnu_lnL'])
        lnL_theta12 = np.array(group['#theta12_lnL @NeutrinoBit::theta12_lnL'])
        lnL_theta13 = np.array(group['#theta13_lnL @NeutrinoBit::theta13_lnL'])
        lnL_theta23 = np.array(group['#theta23_lnL @NeutrinoBit::theta23_lnL'])
        lnL_LUV_LL = np.array(group['#LUV_LL @FlavBit::LUV_likelihood'])
        lnL_sinW2 = np.array(group['#lnL_sinW2 @PrecisionBit::lnL_sinW2_chi2'])

        lnL_list = [lnL_deltaCP, lnL_l2lgamma, lnL_l2lll, lnL_0nubb,
                lnL_W_decays, lnL_W_mass, lnL_Z_inv_width, lnL_bbn,
                lnL_atlas_e, lnL_atlas_mu, lnL_charm_e, lnL_charm_mu,
                lnL_charm_tau, lnL_ckm, lnL_delphi, lnL_e949, lnL_lepuniv,
                lnL_nutev, lnL_pienu, lnL_ps191_e, lnL_ps191_mu, lnL_md21,
                lnL_md3l, lnL_mu2e, lnL_pert, lnL_sum_mnu, lnL_theta12,
                lnL_theta13, lnL_theta23, lnL_LUV_LL, lnL_sinW2]

        lnL_names = ["lnL_deltaCP", "lnL_l2lgamma", "lnL_l2lll", "lnL_0nubb",
               "lnL_W_decays", "lnL_W_mass", "lnL_Z_inv_width", "lnL_bbn",
               "lnL_atlas_e", "lnL_atlas_mu", "lnL_charm_e", "lnL_charm_mu",
               "lnL_charm_tau", "lnL_ckm", "lnL_delphi", "lnL_e949",
               "lnL_lepuniv", "lnL_nutev", "lnL_pienu", "lnL_ps191_e",
               "lnL_ps191_mu", "lnL_md21", "lnL_md3l", "lnL_mu2e", "lnL_pert",
               "lnL_sum_mnu", "lnL_theta12", "lnL_theta13", "lnL_theta23",
               "lnL_LUV_LL", "lnL_sinW2"]

        self.lnL_partial = {}
        for n, l in zip(lnL_names, lnL_list):
            self.lnL_partial[n] = l

        self.lnL = np.array(group['LogLike'])
        self.lnL_sum = sum(lnL_list)

        self.valid = np.array(group['LogLike_isvalid'], dtype = 'bool')

        self.ordering = np.array(group['#ordering @NeutrinoBit::ordering'])

        self.U1 = self.Ue1 + self.Um1 + self.Ut1
        self.U2 = self.Ue2 + self.Um2 + self.Ut2
        self.U3 = self.Ue3 + self.Um3 + self.Ut3
        self.U = self.U1 + self.U2 + self.U3

        self.mNu1 = np.array(group['#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::mNu1'])
        self.mNu2 = np.array(group['#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::mNu2'])
        self.mNu3 = np.array(group['#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::mNu3'])
        self.mMin = np.minimum(np.minimum(self.mNu1, self.mNu2), self.mNu3)

        self.md21 = np.array(group['#md21 @NeutrinoBit::md21'])
        self.md31 = np.array(group['#md31 @NeutrinoBit::md31'])
        self.md32 = np.array(group['#md32 @NeutrinoBit::md32'])

        self.mask2sigma = self.valid & (self.lnL > self.lnL.max() - 0.5*4)

        print "Number of points: %i"%len(self.valid)
        print "Number of valid points: %i"%self.valid.sum()
        print "Number of valid points within 2 sigma: %i"%self.mask2sigma.sum()
        print "Number of valid points within 2 sigma (fraction): %.1e"%(self.mask2sigma.sum()/len(self.valid))

def show_survival_fraction(rhn, sigma = 2):
    for name in rhn.lnL_partial:
        lnL = rhn.lnL_partial[name][rhn.valid]
        lnLmax = lnL.max()
        surv = lnLmax-lnL < 0.5*sigma**2
        print name, surv.sum()/len(lnL)

def show_lnL_hist(rhn, name):
    lnL = rhn.lnL_partial[name][rhn.valid]
    lnL -= lnL.max()
    plt.clf()
    plt.xlabel("renormalized lnL")
    plt.title(name)
    plt.hist(lnL, log=True, bins = 100)
    plt.savefig(OUTPATH + 'hist_'+name+'.pdf')

if __name__ == "__main__":
    rhn = RHN_Chain('/home/cweniger/hdf5/chains/NH_diff.hdf5', 'diff')
    #show_survival_fraction(rhn)
    for name in rhn.lnL_partial:
        print name
        show_lnL_hist(rhn, name)
