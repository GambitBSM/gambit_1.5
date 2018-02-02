#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>

#include "gambit/ColliderBit/analyses/BaseAnalysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"

using namespace std;

/* The CMS 2 lepton direct stop analysis (35.9fb^-1) - `heavy stop'.

   Based on: arXiv:1711.00752
   Yang Zhang

   Known errors:

   Known features:

*/

namespace Gambit {
  namespace ColliderBit {

    class Analysis_CMS_13TeV_2LEPStop_36invfb : public HEPUtilsAnalysis {
    private:

        // Numbers passing cuts
        // int _SRSF[13], _SRDF[13], _SRALL[13],_SRA[3];
        const size_t _SR_size = 13;
        const size_t _SRA_size = 3;
        std::vector<int> _SRSF;
        std::vector<int> _SRDF;
        std::vector<int> _SRALL;
        std::vector<int> _SRA;

        // Cut Flow
        vector<int> cutFlowVector;
        vector<string> cutFlowVector_str;
        int NCUTS;


        // Jet overlap removal
        void JetLeptonOverlapRemoval(vector<HEPUtils::Jet*> &jetvec, vector<HEPUtils::Particle*> &lepvec, double DeltaRMax) {
            //Routine to do jet-lepton check
            //Discards jets if they are within DeltaRMax of a lepton

            vector<HEPUtils::Jet*> Survivors;

            for(unsigned int itjet = 0; itjet < jetvec.size(); itjet++) {
                bool overlap = false;
                HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
                for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
                    HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
                    double dR;

                    dR=jetmom.deltaR_eta(lepmom);

                    if(fabs(dR) <= DeltaRMax) overlap=true;
                }
                if(overlap) continue;
                Survivors.push_back(jetvec.at(itjet));
            }
            jetvec=Survivors;

            return;
        }

        // Lepton overlap removal
        void LeptonJetOverlapRemoval(vector<HEPUtils::Particle*> &lepvec, vector<HEPUtils::Jet*> &jetvec, double DeltaRMax) {
            //Routine to do lepton-jet check
            //Discards leptons if they are within DeltaRMax of a jet

            vector<HEPUtils::Particle*> Survivors;

            for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
                bool overlap = false;
                HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
                for(unsigned int itjet= 0; itjet < jetvec.size(); itjet++) {
                    HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
                    double dR;

                    dR=jetmom.deltaR_eta(lepmom);

                    if(fabs(dR) <= DeltaRMax) overlap=true;
                }
                if(overlap) continue;
                Survivors.push_back(lepvec.at(itlep));
            }
            lepvec=Survivors;

            return;
        }


    public:

        Analysis_CMS_13TeV_2LEPStop_36invfb() {

            for(size_t i=0;i<_SR_size;i++){
                _SRSF[i]=0;
                _SRDF[i]=0;
                _SRALL[i]=0;
            }
            _SRA[0]=0;_SRA[1]=0;_SRA[2]=0;
            NCUTS= 11;
            set_luminosity(35.9);

            for(int i=0;i<NCUTS;i++){
                cutFlowVector.push_back(0);
                cutFlowVector_str.push_back("");
            }

        }

        void analyze(const HEPUtils::Event* event) {
            HEPUtilsAnalysis::analyze(event);

            // Missing energy
            double met = event->met();
            HEPUtils::P4 ptot = event->missingmom();

            // Baseline lepton objects
            const vector<double> a={0,10.};
            const vector<double> b={0,10000.};
            const vector<double> cEl={0.83};
            HEPUtils::BinnedFn2D<double> _eff2dEl(a,b,cEl);
            const vector<double> cMu={0.89};
            HEPUtils::BinnedFn2D<double> _eff2dMu(a,b,cMu);
            vector<HEPUtils::Particle*> baselineElectrons, baselineMuons;
            for (HEPUtils::Particle* electron : event->electrons()) {
                bool hasTrig=has_tag(_eff2dEl, electron->eta(), electron->pT());
                if (electron->pT() > 15. && electron->abseta() < 2.4 && hasTrig) baselineElectrons.push_back(electron);
            }
            for (HEPUtils::Particle* muon : event->muons()) {
                bool hasTrig=has_tag(_eff2dMu, muon->eta(), muon->pT());
                if (muon->pT() > 15. && muon->abseta() < 2.4 && hasTrig) baselineMuons.push_back(muon);
            }
            ATLAS::applyLooseIDElectronSelectionR2(baselineElectrons);
            // Jets
            vector<HEPUtils::Jet*> baselineJets;
            for (HEPUtils::Jet* jet : event->jets()) {
                if (jet->pT() > 30. && fabs(jet->eta()) < 2.4) baselineJets.push_back(jet);
            }

            // Overlap removal
            JetLeptonOverlapRemoval(baselineJets,baselineElectrons,0.4);
            JetLeptonOverlapRemoval(baselineJets,baselineMuons,0.4);

            //Baseline Leptons
            int LooseLepNum = baselineElectrons.size()+baselineMuons.size();
            //Signal Leptons
            ATLAS::applyMediumIDElectronSelectionR2(baselineElectrons);
            vector<HEPUtils::Particle*> signalLeptons;
            for (HEPUtils::Particle* electron : baselineElectrons) {
                signalLeptons.push_back(electron);
            }
            for (HEPUtils::Particle* muon : baselineMuons) {
                signalLeptons.push_back(muon);
            }

            //Put signal jets／leptons in pT order
            //std::sort(signalJets.begin(), signalJets.end(), sortByPT_j);
            //std::sort(signalLeptons.begin(), signalLeptons.end(), sortByPT_l);
            //std::sort(sgJets.begin(), sgJets.end(), sortByPT_j);
            //std::sort(sgLeptons.begin(), sgLeptons.end(), sortByPT_l);

            // Function used to get b jets
            vector<HEPUtils::Jet*> bJets;
            vector<HEPUtils::Jet*> nobJets;
            //const std::vector<double>  a = {0,10.};
            //const std::vector<double>  b = {0,10000.};
            const std::vector<double> c = {0.60};
            HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
            for (HEPUtils::Jet* jet :baselineJets) {
                bool hasTag=has_tag(_eff2d, jet->eta(), jet->pT());
                if(jet->btag() && hasTag && jet->pT() > 25.) {
                        bJets.push_back(jet);
                    }else{
                        nobJets.push_back(jet);
                    }


            }
            int nbjet = bJets.size();
            int njet  = nobJets.size();

            // We now have the signal electrons, muons, jets and b jets- move on to the analysis
            /*********************************************************/
            /*                                                       */
            /* SIGNAL REGIONS                                        */
            /*                                                       */
            /*********************************************************/
            bool cut_2OSLep     =false;
            bool cut_mllGt20    =false;
            bool flg_SF         =false;
            bool cut_mllMZ      =true;
            bool cut_Njet       =false;
            bool cut_Nbjet      =false;
            bool cut_PTmis      =false;
            bool cut_SGt5       =false;
            bool cut_csj1       =false;
            bool cut_csj2       =false;
            bool cut_MT2ll140   =false;
            bool sig_MT2bl_0    =false;
            bool sig_MT2bl_100  =false;
            bool sig_MT2bl_200  =false;
            bool sig_MET_80     =false;
            bool sig_MET_200    =false;
            bool sig_MT2ll_100  =false;
            bool sig_MT2ll_140  =false;
            bool sig_MT2ll_240  =false;
            // Two opposite sign leptons, pT(l1,l2)>25,20GeV
            if(signalLeptons.size() == 2 && LooseLepNum ==2){
                if (signalLeptons[0]->pid()*signalLeptons[1]->pid()<0. && signalLeptons[0]->pT() > 25. && signalLeptons[1]->pT() > 20.){
                    cut_2OSLep = true;
                    /* Calculate variables */
                    // Invariant mass of two leptons
                    HEPUtils::P4 lepton0=signalLeptons.at(0)->mom();
                    HEPUtils::P4 lepton1=signalLeptons.at(1)->mom();
                    double Mll= (lepton0+lepton1).m();
                    // S=MET/sqrt(HT)
                    double HT = 0.;
                    for (HEPUtils::Jet* jet :baselineJets) {
                        HT += jet->pT();
                    }
                    double S=met/sqrt(HT);

                    // Set flags
                    cut_mllGt20 = Mll>20.;
                    flg_SF      = signalLeptons[0]->pid()+signalLeptons[1]->pid()==0;
                    cut_mllMZ   = !(flg_SF && abs(Mll-91.2)<15.);
                    cut_Njet    = njet+nbjet>=2;
                    cut_Nbjet   = nbjet>=1;
                    cut_PTmis   = met>80.;
                    cut_SGt5    = S>5.;

                    // Angular speration of P_T^{miss} and (sub-)leading jet
                    if (cut_Njet) {
                        double cosj1 = cos(baselineJets[0]->phi() - ptot.phi());
                        double cosj2 = cos(baselineJets[1]->phi() - ptot.phi());
                        cut_csj1    = cosj1<0.80;
                        cut_csj2    = cosj2<0.96;
                    }
                    // only calculate mt2 after pass these cuts, to save time
                    if(cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet && cut_PTmis && cut_SGt5 && cut_csj1 && cut_csj2){
                        // MT2
                        double pmiss[3] = { 0, ptot.px(), ptot.py() };
                        mt2_bisect::mt2 mt2_event_bl,mt2_event_ll;
                        // MT2_{ll}
                        double mt2ll=0;
                        double pa_ll[3] = { 0, signalLeptons[0]->mom().px(), signalLeptons[0]->mom().py() };
                        double pb_ll[3] = { 0, signalLeptons[1]->mom().px(), signalLeptons[1]->mom().py() };
                        mt2_event_ll.set_momenta(pa_ll,pb_ll,pmiss);
                        mt2_event_ll.set_mn(0.);
                        mt2ll = mt2_event_ll.get_mt2();
                        // MT2_{blbl}
                        double mt2blbl=0;
                        // find lepton-jet pair minimizes the maximum invariant mass of lepton-jet pairs
                        HEPUtils::P4 bj1 = bJets.at(0)->mom();
                        HEPUtils::P4 bj2;
                        if (nbjet==1) {
                            bj2 = nobJets.at(0)->mom();
                        }else{
                            bj2 = bJets.at(1)->mom();
                        }

                        HEPUtils::P4 l1b1 = lepton0+bj1;
                        HEPUtils::P4 l2b2 = lepton1+bj2;

                        HEPUtils::P4 l1b2 = lepton0+bj2;
                        HEPUtils::P4 l2b1 = lepton1+bj1;
                        double pa_bl[3];
                        double pb_bl[3];
                        pa_bl[0] = 0;
                        pb_bl[0] = 0;
                        if (max(l1b1.m(),l2b2.m()) < max(l1b2.m(),l2b1.m())){
                            pa_bl[1] = l1b1.px();
                            pa_bl[2] = l1b1.py();
                            pb_bl[1] = l2b2.px();
                            pb_bl[2] = l2b2.py();
                        }else{
                            pa_bl[1] = l1b2.px();
                            pa_bl[2] = l1b2.py();
                            pb_bl[1] = l2b1.px();
                            pb_bl[2] = l2b1.py();
                        }
                        mt2_event_bl.set_momenta(pa_bl,pb_bl,pmiss);
                        mt2_event_bl.set_mn(0.);
                        mt2blbl = mt2_event_bl.get_mt2();
                        cut_MT2ll140   = mt2ll>140.;

                        sig_MET_80     = met<200.;
                        sig_MET_200    = met>200.;

                        sig_MT2bl_0    = (mt2blbl<100)&&(mt2blbl>0);
                        sig_MT2bl_100  = (mt2blbl>100)&& (mt2blbl<200);
                        sig_MT2bl_200  = mt2blbl>200;

                        sig_MT2ll_100  = (mt2ll>100.)&&(mt2ll<140.);
                        sig_MT2ll_140  = (mt2ll>140.)&&(mt2ll<240.);
                        sig_MT2ll_240  = (mt2ll>240.);
                    }
                }

            }
            /*********************************************************/
            /*                                                       */
            /* Cut Flow                                              */
            /*                                                       */
            /*********************************************************/
            cutFlowVector_str[0] = "Total ";
            cutFlowVector_str[1] = "2 OS lepton";
            cutFlowVector_str[2] = "m(ll)>20 GeV";
            cutFlowVector_str[3] = "|m(ll)-mZ|>15 GeV";
            cutFlowVector_str[4] = "Njets>2";
            cutFlowVector_str[5] = "Nbjets>1";
            cutFlowVector_str[6] = "MET>80 GeV";
            cutFlowVector_str[7] = "S>5 GeV^{1/2}";
            cutFlowVector_str[8] = "cosPhi(MET,j1)<0.80";
            cutFlowVector_str[9] = "cosPhi(MET,j2)<0.96";
            cutFlowVector_str[10] = "MT2(ll)>140";

            for(int j=0;j<NCUTS;j++){
                if(
                   (j==0) ||
                   (j==1  && cut_2OSLep)||
                   (j==2  && cut_2OSLep && cut_mllGt20)||
                   (j==3  && cut_2OSLep && cut_mllGt20 && cut_mllMZ)||
                   (j==4  && cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet)||
                   (j==5  && cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet)||
                   (j==6  && cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet && cut_PTmis)||
                   (j==7  && cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet && cut_PTmis && cut_SGt5)||
                   (j==8  && cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet && cut_PTmis && cut_SGt5 && cut_csj1) ||
                   (j==9  && cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet && cut_PTmis && cut_SGt5 && cut_csj1 && cut_csj2) ||
                   (j==10 && cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet && cut_PTmis && cut_SGt5 && cut_csj1 && cut_csj2 && cut_MT2ll140)
                   )cutFlowVector[j]++;
            }
            bool pre_cut= cut_2OSLep && cut_mllGt20 && cut_mllMZ && cut_Njet && cut_Nbjet && cut_PTmis && cut_SGt5 && cut_csj1 && cut_csj2 ;
            // signal region          
            for(size_t j=0;j<_SR_size;j++){
                // same flavour
                if(
                   (j==0 && pre_cut && flg_SF && sig_MT2ll_100 && sig_MT2bl_0   && sig_MET_80 )||
                   (j==1 && pre_cut && flg_SF && sig_MT2ll_100 && sig_MT2bl_0   && sig_MET_200)||
                   (j==2 && pre_cut && flg_SF && sig_MT2ll_100 && sig_MT2bl_100 && sig_MET_80 )||
                   (j==3 && pre_cut && flg_SF && sig_MT2ll_100 && sig_MT2bl_100 && sig_MET_200)||
                   (j==4 && pre_cut && flg_SF && sig_MT2ll_100 && sig_MT2bl_200 && sig_MET_80 )||
                   (j==5 && pre_cut && flg_SF && sig_MT2ll_100 && sig_MT2bl_200 && sig_MET_200)||
                   (j==6  && pre_cut && flg_SF && sig_MT2ll_140 && sig_MT2bl_0   && sig_MET_80 )||
                   (j==7  && pre_cut && flg_SF && sig_MT2ll_140 && sig_MT2bl_0   && sig_MET_200)||
                   (j==8  && pre_cut && flg_SF && sig_MT2ll_140 && sig_MT2bl_100 && sig_MET_80 )||
                   (j==9  && pre_cut && flg_SF && sig_MT2ll_140 && sig_MT2bl_100 && sig_MET_200)||
                   (j==10 && pre_cut && flg_SF && sig_MT2ll_140 && sig_MT2bl_200 && sig_MET_80 )||
                   (j==11 && pre_cut && flg_SF && sig_MT2ll_140 && sig_MT2bl_200 && sig_MET_200)||
                   (j==12 && pre_cut && flg_SF && sig_MT2ll_240)
                   )_SRSF[j]++;
                 // diferent flavour
                if(
                   (j==0 && pre_cut && !flg_SF && sig_MT2ll_100 && sig_MT2bl_0   && sig_MET_80 )||
                   (j==1 && pre_cut && !flg_SF && sig_MT2ll_100 && sig_MT2bl_0   && sig_MET_200)||
                   (j==2 && pre_cut && !flg_SF && sig_MT2ll_100 && sig_MT2bl_100 && sig_MET_80 )||
                   (j==3 && pre_cut && !flg_SF && sig_MT2ll_100 && sig_MT2bl_100 && sig_MET_200)||
                   (j==4 && pre_cut && !flg_SF && sig_MT2ll_100 && sig_MT2bl_200 && sig_MET_80 )||
                   (j==5 && pre_cut && !flg_SF && sig_MT2ll_100 && sig_MT2bl_200 && sig_MET_200)||
                   (j==6  && pre_cut && !flg_SF && sig_MT2ll_140 && sig_MT2bl_0   && sig_MET_80 )||
                   (j==7  && pre_cut && !flg_SF && sig_MT2ll_140 && sig_MT2bl_0   && sig_MET_200)||
                   (j==8  && pre_cut && !flg_SF && sig_MT2ll_140 && sig_MT2bl_100 && sig_MET_80 )||
                   (j==9  && pre_cut && !flg_SF && sig_MT2ll_140 && sig_MT2bl_100 && sig_MET_200)||
                   (j==10 && pre_cut && !flg_SF && sig_MT2ll_140 && sig_MT2bl_200 && sig_MET_80 )||
                   (j==11 && pre_cut && !flg_SF && sig_MT2ll_140 && sig_MT2bl_200 && sig_MET_200)||
                   (j==12 && pre_cut && !flg_SF && sig_MT2ll_240)
                   )_SRDF[j]++;
                 // all
                if(
                   (j==0 && pre_cut && sig_MT2ll_100 && sig_MT2bl_0   && sig_MET_80 )||
                   (j==1 && pre_cut && sig_MT2ll_100 && sig_MT2bl_0   && sig_MET_200)||
                   (j==2 && pre_cut && sig_MT2ll_100 && sig_MT2bl_100 && sig_MET_80 )||
                   (j==3 && pre_cut && sig_MT2ll_100 && sig_MT2bl_100 && sig_MET_200)||
                   (j==4 && pre_cut && sig_MT2ll_100 && sig_MT2bl_200 && sig_MET_80 )||
                   (j==5 && pre_cut && sig_MT2ll_100 && sig_MT2bl_200 && sig_MET_200)||
                   (j==6  && pre_cut && sig_MT2ll_140 && sig_MT2bl_0   && sig_MET_80 )||
                   (j==7  && pre_cut && sig_MT2ll_140 && sig_MT2bl_0   && sig_MET_200)||
                   (j==8  && pre_cut && sig_MT2ll_140 && sig_MT2bl_100 && sig_MET_80 )||
                   (j==9  && pre_cut && sig_MT2ll_140 && sig_MT2bl_100 && sig_MET_200)||
                   (j==10 && pre_cut && sig_MT2ll_140 && sig_MT2bl_200 && sig_MET_80 )||
                   (j==11 && pre_cut && sig_MT2ll_140 && sig_MT2bl_200 && sig_MET_200)||
                   (j==12 && pre_cut && sig_MT2ll_240)
                   )_SRALL[j]++;
            }
            for(size_t j=0;j<_SRA_size;j++){
                if(
                   (j==0  && pre_cut && sig_MT2ll_100 && sig_MET_200) ||
                   (j==1  && pre_cut && sig_MT2ll_140 && sig_MET_200)||
                   (j==2  && pre_cut && sig_MT2ll_240)
                   )_SRA[j]++;
            }
        return;

        }


        void add(BaseAnalysis* other) {
            // The base class add function handles the signal region vector and total # events.
            HEPUtilsAnalysis::add(other);

            Analysis_CMS_13TeV_2LEPStop_36invfb* specificOther
                = dynamic_cast<Analysis_CMS_13TeV_2LEPStop_36invfb*>(other);
            // Here we will add the subclass member variables:
            if (NCUTS != specificOther->NCUTS) NCUTS = specificOther->NCUTS;
            for (int j=0; j<NCUTS; j++) {
                cutFlowVector[j] += specificOther->cutFlowVector[j];
                cutFlowVector_str[j] = specificOther->cutFlowVector_str[j];
            }
            for (size_t j=0; j<_SR_size; j++) {
                _SRSF[j] += specificOther->_SRSF[j];
                _SRDF[j] += specificOther->_SRDF[j];
                _SRALL[j] += specificOther->_SRALL[j];
            }
            for (size_t j=0; j<_SRA_size; j++) {
                _SRA[j] += specificOther->_SRA[j];
            }

        }


        void collect_results() {

            double scale_by=1./10000*41.8*35.9;
            cout << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;
            cout << "CUT FLOW: CMS 13 TeV 2 lep stop paper "<<endl;
            cout << "------------------------------------------------------------------------------------------------------------------------------"<<endl;
            cout<< right << setw(40) << "CUT" <<  "," << setw(20) << "RAW" <<  "," << setw(20) << "SCALED"
            <<  "," << setw(20) << "%" <<  "," << setw(20) << "clean adj RAW"<<  "," << setw(20) << "clean adj %" << endl;
            for (int j=0; j<NCUTS; j++) {
                cout << right <<  setw(40) << cutFlowVector_str[j].c_str() <<  "," << setw(20)
                << cutFlowVector[j] <<  "," << setw(20) << cutFlowVector[j]*scale_by <<  "," << setw(20)
                << 100.*cutFlowVector[j]/cutFlowVector[0] << "%" <<  "," << setw(20)
                << cutFlowVector[j]*scale_by <<  "," << setw(20) << 100.*cutFlowVector[j]/cutFlowVector[0]<< "%" << endl;
            }
            for (size_t j=0; j<_SR_size; j++) {
                cout << right <<  setw(40) << "SR_SF_"<<j <<  "," << setw(20)
                << _SRSF[j] <<  "," << setw(20) << _SRSF[j]*scale_by <<  "," << setw(20)
                << 100.*_SRSF[j]/cutFlowVector[0] << "%" <<  "," << setw(20)
                << _SRSF[j]*scale_by <<  "," << setw(20) << 100.*_SRSF[j]/cutFlowVector[0]<< "%" << endl;
            }
           for (size_t j=0; j<_SR_size; j++) {
                cout << right <<  setw(40) << "SR_DF_"<<j <<  "," << setw(20)
                << _SRDF[j] <<  "," << setw(20) << _SRDF[j]*scale_by <<  "," << setw(20)
                << 100.*_SRDF[j]/cutFlowVector[0] << "%" <<  "," << setw(20)
                << _SRDF[j]*scale_by <<  "," << setw(20) << 100.*_SRDF[j]/cutFlowVector[0]<< "%" << endl;
            }
           for (size_t j=0; j<_SR_size; j++) {
                cout << right <<  setw(40) << "SR_ALL_"<<j <<  "," << setw(20)
                << _SRALL[j] <<  "," << setw(20) << _SRALL[j]*scale_by <<  "," << setw(20)
                << 100.*_SRALL[j]/cutFlowVector[0] << "%" <<  "," << setw(20)
                << _SRALL[j]*scale_by <<  "," << setw(20) << 100.*_SRALL[j]/cutFlowVector[0]<< "%" << endl;
            }
           for (size_t j=0; j<_SRA_size; j++) {
                cout << right <<  setw(40) << "SR_A_"<<j <<  "," << setw(20)
                << _SRA[j] <<  "," << setw(20) << _SRA[j]*scale_by <<  "," << setw(20)
                << 100.*_SRA[j]/cutFlowVector[0] << "%" <<  "," << setw(20)
                << _SRA[j]*scale_by <<  "," << setw(20) << 100.*_SRA[j]/cutFlowVector[0]<< "%" << endl;
            }
            cout << "------------------------------------------------------------------------------------------------------------------------------ "<<endl;

            // Same flavour
            /*SignalRegionData results_SRSF0;
            results_SRSF0.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF0.sr_label = "SRSF0";
            results_SRSF0.n_observed = 112.;
            results_SRSF0.n_background = 131.;
            results_SRSF0.background_sys = 30.;
            results_SRSF0.signal_sys = 0.;
            results_SRSF0.n_signal = _SRSF[0];
            add_result(results_SRSF0);

            SignalRegionData results_SRSF1;
            results_SRSF1.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF1.sr_label = "SRSF1";
            results_SRSF1.n_observed = 7.;
            results_SRSF1.n_background = 4.1;
            results_SRSF1.background_sys = 1.1;
            results_SRSF1.signal_sys = 0.;
            results_SRSF1.n_signal = _SRSF[1];
            add_result(results_SRSF1);

            SignalRegionData results_SRSF2;
            results_SRSF2.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF2.sr_label = "SRSF2";
            results_SRSF2.n_observed = 69.;
            results_SRSF2.n_background = 60.;
            results_SRSF2.background_sys = 13.;
            results_SRSF2.signal_sys = 0.;
            results_SRSF2.n_signal = _SRSF[2];
            add_result(results_SRSF2);

            SignalRegionData results_SRSF3;
            results_SRSF3.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF3.sr_label = "SRSF3";
            results_SRSF3.n_observed = 1.;
            results_SRSF3.n_background = 4.8;
            results_SRSF3.background_sys = 1.2;
            results_SRSF3.signal_sys = 0.;
            results_SRSF3.n_signal = _SRSF[3];
            add_result(results_SRSF3);

            SignalRegionData results_SRSF4;
            results_SRSF4.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF4.sr_label = "SRSF4";
            results_SRSF4.n_observed = 0.;
            results_SRSF4.n_background = 0.5;
            results_SRSF4.background_sys = 0.2;
            results_SRSF4.signal_sys = 0.;
            results_SRSF4.n_signal = _SRSF[4];
            add_result(results_SRSF4);

            SignalRegionData results_SRSF5;
            results_SRSF5.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF5.sr_label = "SRSF5";
            results_SRSF5.n_observed = 2.;
            results_SRSF5.n_background = 1.9;
            results_SRSF5.background_sys = 0.5;
            results_SRSF5.signal_sys = 0.;
            results_SRSF5.n_signal = _SRSF[5];
            add_result(results_SRSF5);

            SignalRegionData results_SRSF6;
            results_SRSF6.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF6.sr_label = "SRSF6";
            results_SRSF6.n_observed = 2.;
            results_SRSF6.n_background = 1.1;
            results_SRSF6.background_sys = 0.6;
            results_SRSF6.signal_sys = 0.;
            results_SRSF6.n_signal = _SRSF[6];
            add_result(results_SRSF6);

            SignalRegionData results_SRSF7;
            results_SRSF7.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF7.sr_label = "SRSF7";
            results_SRSF7.n_observed = 2.;
            results_SRSF7.n_background = 0.6;
            results_SRSF7.background_sys = 0.3;
            results_SRSF7.signal_sys = 0.;
            results_SRSF7.n_signal = _SRSF[7];
            add_result(results_SRSF7);

            SignalRegionData results_SRSF8;
            results_SRSF8.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF8.sr_label = "SRSF8";
            results_SRSF8.n_observed = 1.;
            results_SRSF8.n_background = 2.1;
            results_SRSF8.background_sys = 0.7;
            results_SRSF8.signal_sys = 0.;
            results_SRSF8.n_signal = _SRSF[8];
            add_result(results_SRSF8);

            SignalRegionData results_SRSF9;
            results_SRSF9.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF9.sr_label = "SRSF9";
            results_SRSF9.n_observed = 1.;
            results_SRSF9.n_background = 1.6;
            results_SRSF9.background_sys = 0.4;
            results_SRSF9.signal_sys = 0.;
            results_SRSF9.n_signal = _SRSF[9];
            add_result(results_SRSF9);

            SignalRegionData results_SRSF10;
            results_SRSF10.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF10.sr_label = "SRSF10";
            results_SRSF10.n_observed = 0.;
            results_SRSF10.n_background = 0.3;
            results_SRSF10.background_sys = 0.1;
            results_SRSF10.signal_sys = 0.;
            results_SRSF10.n_signal = _SRSF[10];
            add_result(results_SRSF10);

            SignalRegionData results_SRSF11;
            results_SRSF11.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF11.sr_label = "SRSF11";
            results_SRSF11.n_observed = 2.;
            results_SRSF11.n_background = 1.7;
            results_SRSF11.background_sys = 0.4;
            results_SRSF11.signal_sys = 0.;
            results_SRSF11.n_signal = _SRSF[11];
            add_result(results_SRSF11);

            SignalRegionData results_SRSF12;
            results_SRSF12.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRSF12.sr_label = "SRSF12";
            results_SRSF12.n_observed = 1.;
            results_SRSF12.n_background = 0.7;
            results_SRSF12.background_sys = 0.3;
            results_SRSF12.signal_sys = 0.;
            results_SRSF12.n_signal = _SRSF[12];
            add_result(results_SRSF12);

            // Different falvor
            SignalRegionData results_SRDF0;
            results_SRDF0.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF0.sr_label = "SRDF0";
            results_SRDF0.n_observed = 141.;
            results_SRDF0.n_background = 139.;
            results_SRDF0.background_sys = 32.;
            results_SRDF0.signal_sys = 0.;
            results_SRDF0.n_signal = _SRDF[0];
            add_result(results_SRDF0);

            SignalRegionData results_SRDF1;
            results_SRDF1.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF1.sr_label = "SRDF1";
            results_SRDF1.n_observed = 6.;
            results_SRDF1.n_background = 4.0;
            results_SRDF1.background_sys = 1.1;
            results_SRDF1.signal_sys = 0.;
            results_SRDF1.n_signal = _SRDF[1];
            add_result(results_SRDF1);

            SignalRegionData results_SRDF2;
            results_SRDF2.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF2.sr_label = "SRDF2";
            results_SRDF2.n_observed = 67.;
            results_SRDF2.n_background = 70.;
            results_SRDF2.background_sys = 17.;
            results_SRDF2.signal_sys = 0.;
            results_SRDF2.n_signal = _SRDF[2];
            add_result(results_SRDF2);

            SignalRegionData results_SRDF3;
            results_SRDF3.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF3.sr_label = "SRDF3";
            results_SRDF3.n_observed = 5.;
            results_SRDF3.n_background = 3.9;
            results_SRDF3.background_sys = 1.0;
            results_SRDF3.signal_sys = 0.;
            results_SRDF3.n_signal = _SRDF[3];
            add_result(results_SRDF3);

            SignalRegionData results_SRDF4;
            results_SRDF4.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF4.sr_label = "SRDF4";
            results_SRDF4.n_observed = 1.;
            results_SRDF4.n_background = 0.7;
            results_SRDF4.background_sys = 0.2;
            results_SRDF4.signal_sys = 0.;
            results_SRDF4.n_signal = _SRDF[4];
            add_result(results_SRDF4);

            SignalRegionData results_SRDF5;
            results_SRDF5.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF5.sr_label = "SRDF5";
            results_SRDF5.n_observed = 1.;
            results_SRDF5.n_background = 2.1;
            results_SRDF5.background_sys = 0.5;
            results_SRDF5.signal_sys = 0.;
            results_SRDF5.n_signal = _SRDF[5];
            add_result(results_SRDF5);

            SignalRegionData results_SRDF6;
            results_SRDF6.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF6.sr_label = "SRDF6";
            results_SRDF6.n_observed = 1.;
            results_SRDF6.n_background = 0.5;
            results_SRDF6.background_sys = 0.2;
            results_SRDF6.signal_sys = 0.;
            results_SRDF6.n_signal = _SRDF[6];
            add_result(results_SRDF6);

            SignalRegionData results_SRDF7;
            results_SRDF7.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF7.sr_label = "SRDF7";
            results_SRDF7.n_observed = 0.;
            results_SRDF7.n_background = 0.3;
            results_SRDF7.background_sys = 0.2;
            results_SRDF7.signal_sys = 0.;
            results_SRDF7.n_signal = _SRDF[7];
            add_result(results_SRDF7);

            SignalRegionData results_SRDF8;
            results_SRDF8.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF8.sr_label = "SRDF8";
            results_SRDF8.n_observed = 1.;
            results_SRDF8.n_background = 0.8;
            results_SRDF8.background_sys = 0.2;
            results_SRDF8.signal_sys = 0.;
            results_SRDF8.n_signal = _SRDF[8];
            add_result(results_SRDF8);

            SignalRegionData results_SRDF9;
            results_SRDF9.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF9.sr_label = "SRDF9";
            results_SRDF9.n_observed = 0.;
            results_SRDF9.n_background = 0.9;
            results_SRDF9.background_sys = 0.3;
            results_SRDF9.signal_sys = 0.;
            results_SRDF9.n_signal = _SRDF[9];
            add_result(results_SRDF9);

            SignalRegionData results_SRDF10;
            results_SRDF10.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF10.sr_label = "SRDF10";
            results_SRDF10.n_observed = 0.;
            results_SRDF10.n_background = 0.1;
            results_SRDF10.background_sys = 0.1;
            results_SRDF10.signal_sys = 0.;
            results_SRDF10.n_signal = _SRDF[10];
            add_result(results_SRDF10);

            SignalRegionData results_SRDF11;
            results_SRDF11.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF11.sr_label = "SRDF11";
            results_SRDF11.n_observed = 1.;
            results_SRDF11.n_background = 1.2;
            results_SRDF11.background_sys = 0.3;
            results_SRDF11.signal_sys = 0.;
            results_SRDF11.n_signal = _SRDF[11];
            add_result(results_SRDF11);

            SignalRegionData results_SRDF12;
            results_SRDF12.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRDF12.sr_label = "SRDF12";
            results_SRDF12.n_observed = 0.;
            results_SRDF12.n_background = 0.5;
            results_SRDF12.background_sys = 0.2;
            results_SRDF12.signal_sys = 0.;
            results_SRDF12.n_signal = _SRDF[12];
            add_result(results_SRDF12);

            // DF+SF
            SignalRegionData results_SRALL0;
            results_SRALL0.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL0.sr_label = "SRALL0";
            results_SRALL0.n_observed = 253.;
            results_SRALL0.n_background = 271.;
            results_SRALL0.background_sys = 61.;
            results_SRALL0.signal_sys = 0.;
            results_SRALL0.n_signal = _SRALL[0];
            add_result(results_SRALL0);

            SignalRegionData results_SRALL1;
            results_SRALL1.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL1.sr_label = "SRALL1";
            results_SRALL1.n_observed = 13.;
            results_SRALL1.n_background = 8.1;
            results_SRALL1.background_sys = 2.0;
            results_SRALL1.signal_sys = 0.;
            results_SRALL1.n_signal = _SRALL[1];
            add_result(results_SRALL1);

            SignalRegionData results_SRALL2;
            results_SRALL2.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL2.sr_label = "SRALL2";
            results_SRALL2.n_observed = 136.;
            results_SRALL2.n_background = 130.;
            results_SRALL2.background_sys = 29.;
            results_SRALL2.signal_sys = 0.;
            results_SRALL2.n_signal = _SRALL[2];
            add_result(results_SRALL2);

            SignalRegionData results_SRALL3;
            results_SRALL3.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL3.sr_label = "SRALL3";
            results_SRALL3.n_observed = 6.;
            results_SRALL3.n_background = 8.7;
            results_SRALL3.background_sys = 2.0;
            results_SRALL3.signal_sys = 0.;
            results_SRALL3.n_signal = _SRALL[3];
            add_result(results_SRALL3);

            SignalRegionData results_SRALL4;
            results_SRALL4.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL4.sr_label = "SRALL4";
            results_SRALL4.n_observed = 1.;
            results_SRALL4.n_background = 1.2;
            results_SRALL4.background_sys = 0.4;
            results_SRALL4.signal_sys = 0.;
            results_SRALL4.n_signal = _SRALL[4];
            add_result(results_SRALL4);

            SignalRegionData results_SRALL5;
            results_SRALL5.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL5.sr_label = "SRALL5";
            results_SRALL5.n_observed = 3.;
            results_SRALL5.n_background = 4.0;
            results_SRALL5.background_sys = 0.8;
            results_SRALL5.signal_sys = 0.;
            results_SRALL5.n_signal = _SRALL[5];
            add_result(results_SRALL5);

            SignalRegionData results_SRALL6;
            results_SRALL6.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL6.sr_label = "SRALL6";
            results_SRALL6.n_observed = 3.;
            results_SRALL6.n_background = 1.5;
            results_SRALL6.background_sys = 0.7;
            results_SRALL6.signal_sys = 0.;
            results_SRALL6.n_signal = _SRALL[6];
            add_result(results_SRALL6);

            SignalRegionData results_SRALL7;
            results_SRALL7.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL7.sr_label = "SRALL7";
            results_SRALL7.n_observed = 2.;
            results_SRALL7.n_background = 0.8;
            results_SRALL7.background_sys = 0.3;
            results_SRALL7.signal_sys = 0.;
            results_SRALL7.n_signal = _SRALL[7];
            add_result(results_SRALL7);

            SignalRegionData results_SRALL8;
            results_SRALL8.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL8.sr_label = "SRALL8";
            results_SRALL8.n_observed = 2.;
            results_SRALL8.n_background = 2.9;
            results_SRALL8.background_sys = 0.7;
            results_SRALL8.signal_sys = 0.;
            results_SRALL8.n_signal = _SRALL[8];
            add_result(results_SRALL8);

            SignalRegionData results_SRALL9;
            results_SRALL9.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL9.sr_label = "SRALL9";
            results_SRALL9.n_observed = 1.;
            results_SRALL9.n_background = 2.5;
            results_SRALL9.background_sys = 0.5;
            results_SRALL9.signal_sys = 0.;
            results_SRALL9.n_signal = _SRALL[9];
            add_result(results_SRALL9);

            SignalRegionData results_SRALL10;
            results_SRALL10.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL10.sr_label = "SRALL10";
            results_SRALL10.n_observed = 0.;
            results_SRALL10.n_background = 0.4;
            results_SRALL10.background_sys = 0.2;
            results_SRALL10.signal_sys = 0.;
            results_SRALL10.n_signal = _SRALL[10];
            add_result(results_SRALL10);

            SignalRegionData results_SRALL11;
            results_SRALL11.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL11.sr_label = "SRALL11";
            results_SRALL11.n_observed = 3.;
            results_SRALL11.n_background = 2.9;
            results_SRALL11.background_sys = 0.6;
            results_SRALL11.signal_sys = 0.;
            results_SRALL11.n_signal = _SRALL[11];
            add_result(results_SRALL11);

            SignalRegionData results_SRALL12;
            results_SRALL12.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRALL12.sr_label = "SRALL12";
            results_SRALL12.n_observed = 1.;
            results_SRALL12.n_background = 1.1;
            results_SRALL12.background_sys = 0.4;
            results_SRALL12.signal_sys = 0.;
            results_SRALL12.n_signal = _SRALL[12];
            add_result(results_SRALL12);*/

            //Signal region A
            SignalRegionData results_SRA0;
            results_SRA0.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRA0.sr_label = "SRA0";
            results_SRA0.n_observed = 22.;
            results_SRA0.n_background = 20.8;
            results_SRA0.background_sys = 4.4;
            results_SRA0.signal_sys = 0.;
            results_SRA0.n_signal = _SRA[0];
            add_result(results_SRA0);

            SignalRegionData results_SRA1;
            results_SRA1.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRA1.sr_label = "SRA1";
            results_SRA1.n_observed = 6.;
            results_SRA1.n_background = 6.2;
            results_SRA1.background_sys = 1.0;
            results_SRA1.signal_sys = 0.;
            results_SRA1.n_signal = _SRA[1];
            add_result(results_SRA1);

            SignalRegionData results_SRA2;
            results_SRA2.analysis_name = "Analysis_CMS_13TeV_2LEPStop_36invfb";
            results_SRA2.sr_label = "SRA2";
            results_SRA2.n_observed = 1.;
            results_SRA2.n_background = 1.1;
            results_SRA2.background_sys = 0.4;
            results_SRA2.signal_sys = 0.;
            results_SRA2.n_signal = _SRA[2];
            add_result(results_SRA2);

            return;
        }


    protected:
      void clear() {

        std::fill(_SRSF.begin(), _SRSF.end(), 0);
        std::fill(_SRDF.begin(), _SRDF.end(), 0);
        std::fill(_SRALL.begin(), _SRALL.end(), 0);
        std::fill(_SRA.begin(), _SRA.end(), 0);

        std::fill(cutFlowVector.begin(), cutFlowVector.end(), 0);

      }

    };


    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2LEPStop_36invfb)


  }
}
