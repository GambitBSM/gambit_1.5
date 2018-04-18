//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for MicrOmegas SingletDM backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
/// \author Jonathan Cornell
/// \date May 2015, April 2017
///
///  *********************************************

#define BACKENDNAME MicrOmegas_SingletDM
#define BACKENDLANG CC
#define VERSION 3.6.9.2
#define SAFE_VERSION 3_6_9_2

LOAD_LIBRARY

BE_ALLOW_MODELS(SingletDM,SingletDM_running)

BE_FUNCTION(assignVal, int, (char*,double),"assignVal","assignVal", (SingletDM,SingletDM_running))
BE_FUNCTION(vSigma, double, (double, double, int), "vSigma","vSigma", (SingletDM,SingletDM_running))
BE_FUNCTION(darkOmega, double, (double*, int, double), "darkOmega", "oh2", (SingletDM,SingletDM_running))
BE_FUNCTION(sortOddParticles, int, (char*), "sortOddParticles","mass_spectrum", (SingletDM,SingletDM_running))
BE_FUNCTION(cleanDecayTable, void, (), "cleanDecayTable", "cleanDecayTable", (SingletDM,SingletDM_running))
BE_FUNCTION(nucleonAmplitudes, int, (double(*)(double,double,double,double), double*, double*, double*, double*), "nucleonAmplitudes", "nucleonAmplitudes", (SingletDM,SingletDM_running) )
BE_FUNCTION(FeScLoop, double, (double, double, double, double), "FeScLoop", "FeScLoop", (SingletDM,SingletDM_running))
BE_FUNCTION(calcScalarQuarkFF, void, (double, double, double, double), "calcScalarQuarkFF", "calcScalarQuarkFF", (SingletDM,SingletDM_running))

BE_FUNCTION(printChannels, double, (double, double, double, int, FILE*), "printChannels", "momegas_print_channels", (SingletDM,SingletDM_running))

BE_FUNCTION(oneChannel, double, (double,double,char*,char*,char*,char*), "oneChannel", "get_oneChannel", (SingletDM,SingletDM_running))


BE_FUNCTION(mInterp, int, (double,int,int,double*) , "mInterp", "mInterp")
BE_FUNCTION(zInterp, double, (double,double*) , "zInterp", "zInterp")
BE_FUNCTION(readSpectra, int, (), "readSpectra", "readSpectra")

BE_VARIABLE(mocommon_, MicrOmegas::MOcommonSTR, "mocommon_", "MOcommon", (SingletDM,SingletDM_running))
BE_VARIABLE(ForceUG, int, "ForceUG", "ForceUG", (SingletDM,SingletDM_running))
BE_VARIABLE(VZdecay, int, "VZdecay", "VZdecay", (SingletDM,SingletDM_running))
BE_VARIABLE(VWdecay, int, "VWdecay", "VWdecay", (SingletDM,SingletDM_running))

BE_CONV_FUNCTION(dNdE, double, (double,double,int,int), "dNdE")

BE_INI_DEPENDENCY(SingletDM_spectrum, Spectrum)
BE_INI_DEPENDENCY(decay_rates, DecayTable)

#include "gambit/Backends/backend_undefs.hpp"

