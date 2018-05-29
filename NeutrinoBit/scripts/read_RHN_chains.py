#!/usr/bin/env python

from __future__ import division
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

class RHN_Chain(object):
    def __init__(self, INFILE, MODEL, print_keys = False):
        print "Reading %s..."%INFILE
        root = h5py.File(INFILE)
        group = root["RHN"]

        if print_keys:
            for key in group.keys():
                print key
            quit()

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

        lnL_tags = [
            '#deltaCP_lnL @NeutrinoBit::deltaCP_lnL',
            '#l2lgamma_lnL @FlavBit::l2lgamma_likelihood',
            '#l2lll_lnL @FlavBit::l2lll_likelihood',
            '#lnL_0nubb @NeutrinoBit::lnL_0nubb',
            '#lnL_W_decays @PrecisionBit::lnL_W_decays_chi2',
            '#lnL_W_mass @PrecisionBit::lnL_W_mass_chi2',
            '#lnL_Z_inv_width @PrecisionBit::lnL_Z_inv_width_chi2',
            '#lnL_bbn @NeutrinoBit::lnL_bbn',
            '#lnLatlase @NeutrinoBit::lnL_atlas_e',
            '#lnLatlasmu @NeutrinoBit::lnL_atlas_mu',
            '#lnLcharme @NeutrinoBit::lnL_charm_e',
            '#lnLcharmmu @NeutrinoBit::lnL_charm_mu',
            '#lnLcharmtau @NeutrinoBit::lnL_charm_tau',
            '#lnLckm_Vusmin @NeutrinoBit::lnL_ckm_Vusmin',
            '#lnLdelphi @NeutrinoBit::lnL_delphi',
            '#lnLdelphi_shortlived @NeutrinoBit::lnL_delphi_short_lived',
            '#lnLdelphi_longlived @NeutrinoBit::lnL_delphi_long_lived',
            '#lnLe949 @NeutrinoBit::lnL_e949',
            '#lnLlepuniv @NeutrinoBit::lnL_lepuniv',
            '#lnLnutev @NeutrinoBit::lnL_nutev',
            '#lnLpienu @NeutrinoBit::lnL_pienu',
            '#lnLps191e @NeutrinoBit::lnL_ps191_e',
            '#lnLps191mu @NeutrinoBit::lnL_ps191_mu',
            '#md21_lnL @NeutrinoBit::md21_lnL',
            '#md3l_lnL @NeutrinoBit::md3l_lnL',
            '#mu2e_lnL @FlavBit::mu2e_likelihood',
            '#perturbativity_lnL @NeutrinoBit::perturbativity_likelihood',
            '#sum_mnu_lnL @NeutrinoBit::sum_mnu_lnL',
            '#theta12_lnL @NeutrinoBit::theta12_lnL',
            '#theta13_lnL @NeutrinoBit::theta13_lnL',
            '#theta23_lnL @NeutrinoBit::theta23_lnL',
            '#LUV_LL @FlavBit::LUV_likelihood',
            '#lnL_sinW2 @PrecisionBit::lnL_sinW2_chi2',
            '#lnLlhce @NeutrinoBit::lnL_lhc_e',
            '#lnLlhcmu @NeutrinoBit::lnL_lhc_mu']

        lnL_names = ["lnL_deltaCP", "lnL_l2lgamma", "lnL_l2lll", "lnL_0nubb",
               "lnL_W_decays", "lnL_W_mass", "lnL_Z_inv_width", "lnL_bbn",
               "lnL_atlas_e", "lnL_atlas_mu", "lnL_charm_e", "lnL_charm_mu",
               "lnL_charm_tau", "lnL_ckm", "lnL_delphi", "lnL_delphi_short", "lnL_delphi_long", "lnL_e949",
               "lnL_lepuniv", "lnL_nutev", "lnL_pienu", "lnL_ps191_e",
               "lnL_ps191_mu", "lnL_md21", "lnL_md3l", "lnL_mu2e", "lnL_pert",
               "lnL_sum_mnu", "lnL_theta12", "lnL_theta13", "lnL_theta23",
               "lnL_LUV_LL", "lnL_sinW2", "lnL_lhce", "lnL_lhcmu"]

        self.lnL_partial = {}
        for tag, name in zip(lnL_tags, lnL_names):
            try:
                self.lnL_partial[name] = np.array(group[tag])
            except KeyError:
                print "WARNING: Did not find", name
                pass

        self.lnL = np.array(group['LogLike'])

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

