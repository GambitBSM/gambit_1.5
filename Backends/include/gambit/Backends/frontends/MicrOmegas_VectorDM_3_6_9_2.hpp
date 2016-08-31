//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for MicrOmegas VectorDM backend
///  (based on MicrOmegas SingletDM backend)
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
/// \author Ankit Beniwal
/// \date Aug 2016
///
///  *********************************************

#define BACKENDNAME MicrOmegas_VectorDM
#define VERSION 3.6.9.2
#define SAFE_VERSION 3_6_9_2

LOAD_LIBRARY

BE_ALLOW_MODELS(VectorDM, StandardModel_Higgs)

BE_FUNCTION(assignVal, int, (char*,double),"assignVal","assignVal")
BE_FUNCTION(vSigma, double, (double, double, int), "vSigma","vSigma")
BE_FUNCTION(darkOmega, double, (double*, int, double), "darkOmega", "oh2")
BE_FUNCTION(sortOddParticles, int, (char*), "sortOddParticles","mass_spectrum")
BE_FUNCTION(cleanDecayTable, void, (), "cleanDecayTable", "cleanDecayTable")
BE_FUNCTION(calcSpectrum, double, (int, double*, double*, double*, double*, double*, double*, int*), "calcSpectrum", "calcSpectrum")
BE_FUNCTION(nucleonAmplitudes, int, (double(*)(double,double,double,double), double*, double*, double*, double*), "nucleonAmplitudes", "nucleonAmplitudes" )
BE_FUNCTION(FeScLoop, double, (double, double, double, double), "FeScLoop", "FeScLoop")
BE_FUNCTION(calcScalarQuarkFF, void, (double, double, double, double), "calcScalarQuarkFF", "calcScalarQuarkFF")

BE_FUNCTION(mInterp, int, (double,int,int,double*) , "mInterp", "mInterp")
BE_FUNCTION(zInterp, double, (double,double*) , "zInterp", "zInterp")
BE_FUNCTION(readSpectra, int, (), "readSpectra", "readSpectra")

BE_VARIABLE(mocommon_, MicrOmegas::MOcommonSTR, "mocommon_", "MOcommon")
BE_VARIABLE(vSigmaCh, MicrOmegas::aChannel*, "vSigmaCh", "vSigmaCh")
BE_VARIABLE(ForceUG, int, "ForceUG", "ForceUG")
BE_VARIABLE(VZdecay, int, "VZdecay", "VZdecay")
BE_VARIABLE(VWdecay, int, "VWdecay", "VWdecay")

namespace Gambit
{
  namespace Backends
  {
    namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
    {
      double dNdE(double Ecm, double E, int inP, int outN)
      {
        // outN 0-5: gamma, e+, p-, nu_e, nu_mu, nu_tau
        // inP:  0 - 6: glu, d, u, s, c, b, t
        //       7 - 9: e, m, l
        //       10 - 15: Z, ZT, ZL, W, WT, WL
        double tab[250];  // NZ = 250
        readSpectra();
        mInterp(Ecm/2, inP, outN, tab);
        return zInterp(log(E/Ecm*2), tab);
      }

    } /* end namespace BACKENDNAME_SAFE_VERSION */
  } /* end namespace Backends */
} /* end namespace Gambit */

BE_CONV_FUNCTION(dNdE, double, (double,double,int,int), "dNdE")

BE_INI_DEPENDENCY(SMINPUTS, SMInputs)

BE_INI_FUNCTION
{
     int error;
     char cdmName[10];

     // Currently only works correctly in unitary gauge
     *ForceUG=1;

     // Set VectorDM model parameters in micrOmegas
     error = assignVal((char*)"mV", *Param["mV"]);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mV in"
             "MicrOmegas. MicrOmegas error code: " + std::to_string(error));     

     error = assignVal((char*)"lhV", *Param["lambda_hV"]);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set lambda_hV in"
             "MicrOmegas. MicrOmegas error code: " + std::to_string(error));     

     // Set SM + Higgs mass parameters in micrOmegas
     const SMInputs& sminputs = *Dep::SMINPUTS;	

     // EE = sqrt(4*pi*(1/alphainv))
     error = assignVal((char*)"EE", sqrt(4*M_PI*1/(sminputs.alphainv)));
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set alphainv in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // GG = sqrt(4*pi*alphaS)
     // error = assignVal((char*)"GG", sqrt(4*M_PI*sminputs.alphaS));
     //     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set GG in"
     //      " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // v0 = 1/sqrt(sqrt(2)*GF)
     error = assignVal((char*)"v0", 1/sqrt(sqrt(2)*sminputs.GF));
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set GF in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // mu(2 GeV) in MSbar scheme
     error = assignVal((char*)"Mu", sminputs.mU);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mU in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // md(2 GeV) in MSbar scheme
     error = assignVal((char*)"Md", sminputs.mD);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mD in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // ms(2 GeV) in MSbar scheme
     error = assignVal((char*)"Ms", sminputs.mS);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mS in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // mc(mc) in MSbar scheme
     error = assignVal((char*)"Mc", sminputs.mCmC);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mCmC in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // mb(mb) in MSbar scheme
     error = assignVal((char*)"Mb", sminputs.mBmB);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mBmB in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // mtop(pole)
     error = assignVal((char*)"Mtop", sminputs.mT);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mT in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // mE(pole)
     error = assignVal((char*)"Me", sminputs.mE);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mE in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));
             
     // mMu(pole)
     error = assignVal((char*)"Mm", sminputs.mMu);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mMu in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // mTau(pole)
     error = assignVal((char*)"Mtau", sminputs.mTau);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mTau in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));

     // mZ(pole)
     error = assignVal((char*)"MZ", sminputs.mZ);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mZ in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));
     	
     // mh
     error = assignVal((char*)"MH", *Param["mH"]);
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "Unable to set mH in"
             " MicrOmegas. MicrOmegas error code: " + std::to_string(error));
		
     error = sortOddParticles(byVal(cdmName));
     if (error != 0) BackendIniBit_error().raise(LOCAL_INFO, "MicrOmegas function "
             "sortOddParticles returned error code: " + std::to_string(error));

}
END_BE_INI_FUNCTION

#include "gambit/Backends/backend_undefs.hpp"

