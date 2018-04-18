//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for MicrOmegas SingletDMZ3 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
/// \author Jonathan Cornell
/// \date May 2015, April 2017
///
///  *********************************************

#define BACKENDNAME MicrOmegas_SingletDMZ3
#define BACKENDLANG CC
#define VERSION 3.6.9.2
#define SAFE_VERSION 3_6_9_2

LOAD_LIBRARY

BE_ALLOW_MODELS(SingletDMZ3)

BE_FUNCTION(assignVal, int, (char*,double),"assignVal","assignVal", (SingletDMZ3))
BE_FUNCTION(vSigma, double, (double, double, int), "vSigma","vSigma", (SingletDMZ3))
BE_FUNCTION(darkOmega, double, (double*, int, double), "darkOmega", "oh2", (SingletDMZ3))
BE_FUNCTION(sortOddParticles, int, (char*), "sortOddParticles","mass_spectrum", (SingletDMZ3))
BE_FUNCTION(cleanDecayTable, void, (), "cleanDecayTable", "cleanDecayTable", (SingletDMZ3))
BE_FUNCTION(nucleonAmplitudes, int, (double(*)(double,double,double,double), double*, double*, double*, double*), "nucleonAmplitudes", "nucleonAmplitudes", (SingletDMZ3) )
BE_FUNCTION(FeScLoop, double, (double, double, double, double), "FeScLoop", "FeScLoop", (SingletDMZ3))
BE_FUNCTION(calcScalarQuarkFF, void, (double, double, double, double), "calcScalarQuarkFF", "calcScalarQuarkFF", (SingletDMZ3))

BE_FUNCTION(calcSpectrum, double, (int, double*, double*, double*, double*, double*, double*, int*), "calcSpectrum", "calcSpectrum", (SingletDMZ3))


BE_FUNCTION(printChannels, double, (double, double, double, int, FILE*), "printChannels", "momegas_print_channels", (SingletDMZ3))

BE_FUNCTION(oneChannel, double, (double,double,char*,char*,char*,char*), "oneChannel", "get_oneChannel", (SingletDMZ3))


BE_FUNCTION(mInterp, int, (double,int,int,double*) , "mInterp", "mInterp")
BE_FUNCTION(zInterp, double, (double,double*) , "zInterp", "zInterp")
BE_FUNCTION(readSpectra, int, (), "readSpectra", "readSpectra")


BE_VARIABLE(vSigmaCh_, MicrOmegas::aChannel, "_vSigmaCh", "vSigmaCh", (SingletDMZ3))
BE_VARIABLE(mocommon_, MicrOmegas::MOcommonSTR, "mocommon_", "MOcommon", (SingletDMZ3))
BE_VARIABLE(ForceUG, int, "ForceUG", "ForceUG", (SingletDMZ3))
BE_VARIABLE(VZdecay, int, "VZdecay", "VZdecay", (SingletDMZ3))
BE_VARIABLE(VWdecay, int, "VWdecay", "VWdecay", (SingletDMZ3))

BE_CONV_FUNCTION(dNdE, double, (double,double,int,int), "dNdE")

BE_INI_DEPENDENCY(SingletDMZ3_spectrum, Spectrum)
BE_INI_DEPENDENCY(decay_rates, DecayTable)

#include "gambit/Backends/backend_undefs.hpp"

