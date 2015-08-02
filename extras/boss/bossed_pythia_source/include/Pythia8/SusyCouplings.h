#ifndef __boss__SusyCouplings_Pythia_8_209_h__
#define __boss__SusyCouplings_Pythia_8_209_h__

// SusyCouplings.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for setup of common SUSY couplings.

#ifndef Pythia8_SusyCouplings_H
#define Pythia8_SusyCouplings_H

#include "Pythia8/PythiaComplex.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyLesHouches.h"

namespace Pythia8 {

class ParticleData;

//==========================================================================

// CoupSUSY
// Auxiliary class to compute and store various SM and SUSY couplings.

} 
#define ENUMS_DECLARED
#include "backend_types/Pythia_8_209/abstract_CoupSUSY.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
namespace Pythia8 { 
class CoupSUSY : public virtual Abstract_CoupSUSY, public Couplings{

public:

  // Constructor
  CoupSUSY() {isInit=false; isNMSSM = false; isSUSY=true;}

  // Initialize
  void initSUSY(SusyLesHouches* slhaPtrIn, Info* infoPtrIn,
                ParticleData* particleDataPtrIn, Settings* settingsPtrIn);

  // Status flag. Flag for NMSSM.
  bool isInit, isNMSSM;

  // Z and W pole masses and widths
  double mWpole, wWpole, mZpole, wZpole;

  // Running masses and weak mixing angle
  // (default to pole values if no running available)
  double mW, mZ, sin2W, sinW, cosW;

  // Tanbeta
  double tanb, cosb, sinb;

  //Higgs-sector parameters
  double muHiggs, alphaHiggs, mAHiggs;

  // ~qq~g couplings
  complex LsddG[7][4], RsddG[7][4];
  complex LsuuG[7][4], RsuuG[7][4];
  // Assume generation index for Squark. Translate if PDG code instead.
  complex getLsqqG(int iGenSq, int idQ) {if (abs(iGenSq) > 1000000)
      iGenSq =  3*(abs(iGenSq)/2000000) + (abs(iGenSq)%10+1)/2;
    return (abs(idQ)%2 == 0) ? LsuuG[iGenSq][abs(idQ)/2]
      : LsddG[iGenSq][(abs(idQ)+1)/2] ;}
  complex getRsqqG(int iGenSq, int idQ) {if (abs(iGenSq) > 1000000)
      iGenSq =  3*(abs(iGenSq)/2000000) + (abs(iGenSq)%10+1)/2;
    return (abs(idQ)%2 == 0) ? RsuuG[iGenSq][abs(idQ)/2]
      : RsddG[iGenSq][(abs(idQ)+1)/2] ;}

  // ~chi0~chi0Z couplings
  complex OLpp[6][6], ORpp[6][6];

  // ~chi+~chi-Z couplings
  complex OLp[3][3], ORp[3][3];

  // ~chi0~chi+W couplings
  complex OL[6][3], OR[6][3];

  // qqZ couplings
  double LqqZ[7], RqqZ[7];

  // ~q~qZ couplings
  complex LsdsdZ[7][7], RsdsdZ[7][7];
  complex LsusuZ[7][7], RsusuZ[7][7];
  complex getLsqsqZ(int idSq1, int idSq2) {
    if (abs(idSq1)%2 != abs(idSq2)%2) return complex(0.0,0.0);
    int iGen1 = 3*(abs(idSq1)/2000000) + (abs(idSq1)%10+1)/2;
    int iGen2 = 3*(abs(idSq2)/2000000) + (abs(idSq2)%10+1)/2;
    return (abs(idSq1)%2 == 0) ? LsusuZ[iGen1][iGen2] : LsdsdZ[iGen1][iGen2];}
  complex getRsqsqZ(int idSq1, int idSq2) {
    if (abs(idSq1)%2 != abs(idSq2)%2) return complex(0.0,0.0);
    int iGen1 = 3*(abs(idSq1)/2000000) + (abs(idSq1)%10+1)/2;
    int iGen2 = 3*(abs(idSq2)/2000000) + (abs(idSq2)%10+1)/2;
    return (abs(idSq1)%2 == 0) ? RsusuZ[iGen1][iGen2] : RsdsdZ[iGen1][iGen2];}

  // udW couplings
  complex LudW[4][4], RudW[4][4];

  // ~u~dW couplings
  complex LsusdW[7][7], RsusdW[7][7];

  // ~qq~chi0 couplings
  complex LsddX[7][4][6], RsddX[7][4][6];
  complex LsuuX[7][4][6], RsuuX[7][4][6];
  complex getLsqqX(int iSq, int idQ, int iNeut) {return (abs(idQ)%2 == 0)
    ? LsuuX[iSq][abs(idQ)/2][iNeut] : LsddX[iSq][(abs(idQ)+1)/2][iNeut] ;}
  complex getRsqqX(int iSq, int idQ, int iNeut) {return (abs(idQ)%2 == 0)
    ? RsuuX[iSq][abs(idQ)/2][iNeut] : RsddX[iSq][(abs(idQ)+1)/2][iNeut] ;}

  // ~du~chi+ couplings
  complex LsduX[7][4][3], RsduX[7][4][3];

  // ~ud~chi+ couplings
  complex LsudX[7][4][3], RsudX[7][4][3];

  //llZ couplings
  double LllZ[7], RllZ[7];

  //lvW couplings
  complex LlvW[4][4], RlvW[4][4];

  // ~l~lZ couplings
  complex LslslZ[7][7],RslslZ[7][7];
  complex LsvsvZ[7][7],RsvsvZ[7][7];

  // ~l~vW couplings
  complex LslsvW[7][7], RslsvW[7][7];

  // ~ll~chi0 couplings
  complex LsvvX[7][4][6], RsvvX[7][4][6];
  complex LsllX[7][4][6], RsllX[7][4][6];

  // ~vl~chi+ couplings
  complex LsvlX[7][4][3], RsvlX[7][4][3];

  // ~lv~chi+ couplings
  complex LslvX[7][4][3], RslvX[7][4][3];

  // RPV couplings
  double rvLLE[4][4][4], rvLQD[4][4][4], rvUDD[4][4][4];
  // Flags for RPV couplings
  bool isLLE, isLQD, isUDD;

  //Squark and slepton mixing matrix: needed for RPV
  complex Rusq[7][7], Rdsq[7][7];
  complex Rsl[7][7], Rsv[7][7];

  // Return neutralino, chargino, sup, sdown and slepton flavour codes.
  int idNeut(int idChi);
  int idChar(int idChi);
  int idSup(int iSup);
  int idSdown(int iSdown);
  int idSlep(int iSlep);

  //Reverse lookup for neutralinos and charginos
  int typeNeut(int idPDG);
  int typeChar(int idPDG);

  // Pointer to SLHA instance
  // Used in SusyResonanceWidths for checking if decay table exists
  SusyLesHouches* slhaPtr;

private:
  // Debug flag
  static const bool DBSUSY;

  // Pointer to the info class
  Info*          infoPtr;

  // Pointer to the settings database.
  Settings*      settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;


        public:
            Abstract_CoupSUSY* pointerCopy__BOSS();

            void pointerAssign__BOSS(Abstract_CoupSUSY* in);

        public:
            bool& isInit_ref__BOSS();

            bool& isNMSSM_ref__BOSS();

            double& mWpole_ref__BOSS();

            double& wWpole_ref__BOSS();

            double& mZpole_ref__BOSS();

            double& wZpole_ref__BOSS();

            double& mW_ref__BOSS();

            double& mZ_ref__BOSS();

            double& sin2W_ref__BOSS();

            double& sinW_ref__BOSS();

            double& cosW_ref__BOSS();

            double& tanb_ref__BOSS();

            double& cosb_ref__BOSS();

            double& sinb_ref__BOSS();

            double& muHiggs_ref__BOSS();

            double& alphaHiggs_ref__BOSS();

            double& mAHiggs_ref__BOSS();

            std::complex<double> (&LsddG_ref__BOSS())[7][4];

            std::complex<double> (&RsddG_ref__BOSS())[7][4];

            std::complex<double> (&LsuuG_ref__BOSS())[7][4];

            std::complex<double> (&RsuuG_ref__BOSS())[7][4];

            std::complex<double> (&OLpp_ref__BOSS())[6][6];

            std::complex<double> (&ORpp_ref__BOSS())[6][6];

            std::complex<double> (&OLp_ref__BOSS())[3][3];

            std::complex<double> (&ORp_ref__BOSS())[3][3];

            std::complex<double> (&OL_ref__BOSS())[6][3];

            std::complex<double> (&OR_ref__BOSS())[6][3];

            double (&LqqZ_ref__BOSS())[7];

            double (&RqqZ_ref__BOSS())[7];

            std::complex<double> (&LsdsdZ_ref__BOSS())[7][7];

            std::complex<double> (&RsdsdZ_ref__BOSS())[7][7];

            std::complex<double> (&LsusuZ_ref__BOSS())[7][7];

            std::complex<double> (&RsusuZ_ref__BOSS())[7][7];

            std::complex<double> (&LudW_ref__BOSS())[4][4];

            std::complex<double> (&RudW_ref__BOSS())[4][4];

            std::complex<double> (&LsusdW_ref__BOSS())[7][7];

            std::complex<double> (&RsusdW_ref__BOSS())[7][7];

            std::complex<double> (&LsddX_ref__BOSS())[7][4][6];

            std::complex<double> (&RsddX_ref__BOSS())[7][4][6];

            std::complex<double> (&LsuuX_ref__BOSS())[7][4][6];

            std::complex<double> (&RsuuX_ref__BOSS())[7][4][6];

            std::complex<double> (&LsduX_ref__BOSS())[7][4][3];

            std::complex<double> (&RsduX_ref__BOSS())[7][4][3];

            std::complex<double> (&LsudX_ref__BOSS())[7][4][3];

            std::complex<double> (&RsudX_ref__BOSS())[7][4][3];

            double (&LllZ_ref__BOSS())[7];

            double (&RllZ_ref__BOSS())[7];

            std::complex<double> (&LlvW_ref__BOSS())[4][4];

            std::complex<double> (&RlvW_ref__BOSS())[4][4];

            std::complex<double> (&LslslZ_ref__BOSS())[7][7];

            std::complex<double> (&RslslZ_ref__BOSS())[7][7];

            std::complex<double> (&LsvsvZ_ref__BOSS())[7][7];

            std::complex<double> (&RsvsvZ_ref__BOSS())[7][7];

            std::complex<double> (&LslsvW_ref__BOSS())[7][7];

            std::complex<double> (&RslsvW_ref__BOSS())[7][7];

            std::complex<double> (&LsvvX_ref__BOSS())[7][4][6];

            std::complex<double> (&RsvvX_ref__BOSS())[7][4][6];

            std::complex<double> (&LsllX_ref__BOSS())[7][4][6];

            std::complex<double> (&RsllX_ref__BOSS())[7][4][6];

            std::complex<double> (&LsvlX_ref__BOSS())[7][4][3];

            std::complex<double> (&RsvlX_ref__BOSS())[7][4][3];

            std::complex<double> (&LslvX_ref__BOSS())[7][4][3];

            std::complex<double> (&RslvX_ref__BOSS())[7][4][3];

            double (&rvLLE_ref__BOSS())[4][4][4];

            double (&rvLQD_ref__BOSS())[4][4][4];

            double (&rvUDD_ref__BOSS())[4][4][4];

            bool& isLLE_ref__BOSS();

            bool& isLQD_ref__BOSS();

            bool& isUDD_ref__BOSS();

            std::complex<double> (&Rusq_ref__BOSS())[7][7];

            std::complex<double> (&Rdsq_ref__BOSS())[7][7];

            std::complex<double> (&Rsl_ref__BOSS())[7][7];

            std::complex<double> (&Rsv_ref__BOSS())[7][7];



        public:
            void initSUSY__BOSS(Pythia8::Abstract_SusyLesHouches*, Pythia8::Abstract_Info*, Pythia8::Abstract_ParticleData*, Pythia8::Abstract_Settings*);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SusyCouplings_H

#endif /* __boss__SusyCouplings_Pythia_8_209_h__ */
