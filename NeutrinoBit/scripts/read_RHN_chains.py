#!/usr/bin/env python

from __future__ import division
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

class RHN_Chain(object):
    def __init__(self, INFILE, MODEL = 'auto', print_keys = False, renormalize
            = True, sub_slide = True):
        print "Reading %s..."%INFILE
        if MODEL == 'auto':
            if 'diff' in INFILE:
                MODEL = 'diff'
            elif 'full' in INFILE:
                MODEL = 'full'
            else:
                raise KeyError("Model not known.")
        root = h5py.File(INFILE)
        group = root["RHN"]

        if print_keys:
            for key in group.keys():
                print key
            quit()

        SM_mode = 'mNudiff' if any(['mNudiff' in key for key in group.keys()]) else 'SLHA2'

        valid = np.array(group['LogLike_isvalid'], dtype = 'bool')
        def get_data(tag):
            return np.array(group[tag])[valid]

        if MODEL == 'diff':
            self.M1 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_1')
            self.dM21 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::delta_M21')
            self.M3 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::M_3')
            self.ReOmega12 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ReOm12')
            self.ReOmega13 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ReOm13')
            self.ReOmega23 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ReOm23')
            self.ImOmega12 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ImOm12')
            self.ImOmega13 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ImOm13')
            self.ImOmega23 = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::ImOm23')
            self.Rorder = get_data('#RightHandedNeutrinos_diff_parameters @RightHandedNeutrinos_diff::primary_parameters::Rorder')
            self.M2 = self.M1 + self.dM21
        elif MODEL == 'full':
            self.M1 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::M_1')
            self.M2 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::M_2')
            self.M3 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::M_3')
            self.ReOmega12 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ReOm12')
            self.ReOmega13 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ReOm13')
            self.ReOmega23 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ReOm23')
            self.ImOmega12 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ImOm12')
            self.ImOmega13 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ImOm13')
            self.ImOmega23 = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::ImOm23')
            self.Rorder = get_data('#RightHandedNeutrinos_parameters @RightHandedNeutrinos::primary_parameters::Rorder')
            self.dM2 = self.M2 - self.M1
            self.dM3 = self.M3 - self.M1
        else:
            raise KeyError("Model unknown")

        self.Ue1 = get_data('#Ue1 @NeutrinoBit::Ue1')
        self.Ue2 = get_data('#Ue2 @NeutrinoBit::Ue2')
        self.Ue3 = get_data('#Ue3 @NeutrinoBit::Ue3')
        self.Um1 = get_data('#Um1 @NeutrinoBit::Um1')
        self.Um2 = get_data('#Um2 @NeutrinoBit::Um2')
        self.Um3 = get_data('#Um3 @NeutrinoBit::Um3')
        self.Ut1 = get_data('#Ut1 @NeutrinoBit::Ut1')
        self.Ut2 = get_data('#Ut2 @NeutrinoBit::Ut2')
        self.Ut3 = get_data('#Ut3 @NeutrinoBit::Ut3')

        self.Ue1_phase = get_data('#Ue1_phase @NeutrinoBit::Ue1_phase')
        self.Ue2_phase = get_data('#Ue2_phase @NeutrinoBit::Ue2_phase')
        self.Ue3_phase = get_data('#Ue3_phase @NeutrinoBit::Ue3_phase')
        self.Um1_phase = get_data('#Um1_phase @NeutrinoBit::Um1_phase')
        self.Um2_phase = get_data('#Um2_phase @NeutrinoBit::Um2_phase')
        self.Um3_phase = get_data('#Um3_phase @NeutrinoBit::Um3_phase')
        self.Ut1_phase = get_data('#Ut1_phase @NeutrinoBit::Ut1_phase')
        self.Ut2_phase = get_data('#Ut2_phase @NeutrinoBit::Ut2_phase')
        self.Ut3_phase = get_data('#Ut3_phase @NeutrinoBit::Ut3_phase')

        self.theta13 = get_data('#theta13 @NeutrinoBit::theta13')
        self.theta12 = get_data('#theta12 @NeutrinoBit::theta12')
        self.theta23 = get_data('#theta23 @NeutrinoBit::theta23')
        self.deltaCP = get_data('#deltaCP @NeutrinoBit::deltaCP')
        self.alpha1 = get_data('#StandardModel_mNudiff_parameters @StandardModel_mNudiff::primary_parameters::alpha1')
        self.alpha2 = get_data('#StandardModel_mNudiff_parameters @StandardModel_mNudiff::primary_parameters::alpha2')

        self.ue1 = self.Ue1**0.5 * np.exp(1j*self.Ue1_phase)
        self.ue2 = self.Ue2**0.5 * np.exp(1j*self.Ue2_phase)
        self.ue3 = self.Ue3**0.5 * np.exp(1j*self.Ue3_phase)
        self.um1 = self.Um1**0.5 * np.exp(1j*self.Um1_phase)
        self.um2 = self.Um2**0.5 * np.exp(1j*self.Um2_phase)
        self.um3 = self.Um3**0.5 * np.exp(1j*self.Um3_phase)
        self.ut1 = self.Ut1**0.5 * np.exp(1j*self.Ut1_phase)
        self.ut2 = self.Ut2**0.5 * np.exp(1j*self.Ut2_phase)
        self.ut3 = self.Ut3**0.5 * np.exp(1j*self.Ut3_phase)

        try:
            self.mbb = get_data('#mbb_0nubb_Xe @NeutrinoBit::RHN_mbb_0nubb_Xe')
        except KeyError:
            pass

        lnL_tags = [
            '#deltaCP_lnL @NeutrinoBit::deltaCP_lnL',
            '#l2lgamma_lnL @FlavBit::l2lgamma_likelihood',
            '#l2lll_lnL @FlavBit::l2lll_likelihood',
            '#lnL_0nubb @NeutrinoBit::lnL_0nubb',
            '#lnL_mbb_0nubb_KamLAND_Zen @NeutrinoBit::lnL_mbb_0nubb_KamLAND_Zen',
            '#lnL_mbb_0nubb_GERDA @NeutrinoBit::lnL_mbb_0nubb_GERDA',
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

        lnL_names = ["lnL_deltaCP", "lnL_l2lgamma", "lnL_l2lll", 
                "lnL_0nubb", "lnL_mbb_0nubb_KamLAND_Zen", "lnL_mbb_0nubb_GERDA",
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
                self.lnL_partial[name] = get_data(tag)
            except KeyError:
                print "WARNING: Did not find", name
                pass
        
        try:
            self.lnL_slide = get_data("#RHN_coupling_slide @NeutrinoBit::coupling_slide")
            print "Found coupling slide!"
            print self.lnL_slide/np.log10(self.Ue1)
            print self.lnL_slide/np.log10(self.Um1)
            print self.lnL_slide/np.log10(self.Ut1)
            print self.lnL_slide/np.log10(self.Ue1*self.Um1*self.Ut1)
        except KeyError:
            self.lnL_slide = None

        self.lnL = get_data('LogLike')
        if self.lnL_slide is not None and sub_slide:
            self.lnL -= self.lnL_slide

        try:
            self.Vus = get_data('#calc_Vus @NeutrinoBit::calc_Vus')
        except KeyError:
            pass

        #self.ordering = get_data('#ordering @NeutrinoBit::ordering')

        self.U1 = self.Ue1 + self.Um1 + self.Ut1
        self.U2 = self.Ue2 + self.Um2 + self.Ut2
        self.U3 = self.Ue3 + self.Um3 + self.Ut3
        self.U = self.U1 + self.U2 + self.U3

        if SM_mode == 'SLHA2':
            self.mNu1 = get_data('#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::mNu1')
            self.mNu2 = get_data('#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::mNu2')
            self.mNu3 = get_data('#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::mNu3')
        elif SM_mode == 'mNudiff':
            mNu_light = get_data('#StandardModel_mNudiff_parameters @StandardModel_mNudiff::primary_parameters::mNu_light')
            dmNu21 = get_data('#StandardModel_mNudiff_parameters @StandardModel_mNudiff::primary_parameters::dmNu21')
            dmNu3l = get_data('#StandardModel_mNudiff_parameters @StandardModel_mNudiff::primary_parameters::dmNu3l')
            self.mNu1 = np.where(dmNu3l > 0, mNu_light, (mNu_light**2 + abs(dmNu3l))**0.5)*1e-9
            self.mNu2 = np.where(dmNu3l > 0, (mNu_light**2 + dmNu21)**0.5, (mNu_light**2 + abs(dmNu3l) + dmNu21)**0.5)*1e-9
            self.mNu3 = np.where(dmNu3l > 0, (mNu_light**2 + dmNu21 + abs(dmNu3l))**0.5, mNu_light)*1e-9
        else:
            raise KeyError()
            

        self.mMin = np.minimum(np.minimum(self.mNu1, self.mNu2), self.mNu3)
        self.md21 = get_data('#md21 @NeutrinoBit::md21')
        self.md31 = get_data('#md31 @NeutrinoBit::md31')
        self.md32 = get_data('#md32 @NeutrinoBit::md32')

        if renormalize:
            self.lnL[self.md31>0] -= self.lnL[self.md31>0].max()
            self.lnL[self.md31<0] -= self.lnL[self.md31<0].max()

        mask2sigma = self.lnL > self.lnL.max() - 0.5*4

        print "Number of points: %i"%len(valid)
        print "Number of valid points: %i"%valid.sum()
        print "Number of valid points within 2 sigma: %i"%mask2sigma.sum()
        print "Number of valid points within 2 sigma (fraction): %.1e"%(mask2sigma.sum()/len(valid))

