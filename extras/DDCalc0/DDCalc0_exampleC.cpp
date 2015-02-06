/**********************************************************************
 * DDCALC0 EXAMPLE PROGRAM (C++)
 * This program shows how to use the DDCalc0 module from C++, making
 * use of the interface defined in the DDCalc0.hpp header file.
 * 
 * Compile (assuming DDCalc0.o is already compiled):
 *   g++ -o DDCalc0_exampleC DDCalc0_exampleC.cpp DDCalc0.o -lgfortran
 *   icc -o DDCalc0_exampleC DDCalc0_exampleC.cpp DDCalc0.o -lifcore
 *
 * Run:
 *   ./DDCalc0_exampleC [--mG|--mfa|--msigma]
 * where the optional flag specifies the form in which the WIMP-nucleon
 * couplings will be provided (default: --msigma).
 * 
 * 
 *       A. Scaffidi     U of Adelaide    2015    
 *       C. Savage       Nordita          2015
 * 
 **********************************************************************/


/* All the DDCalc0 routines used below are declared in the DDCalc0.hpp
   header file. */

////////////////////////
#include "DDCalc0.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;
////////////////////////


// CONSTANTS -----------------------------------------------------------

// These constants will be used to specify the type of input parameters.
const int TYPE_MG     = 1;  // Four-fermion effective couplings G
const int TYPE_MFA    = 2;  // Effective couplings f (SI), a (SD)
const int TYPE_MSIGMA = 3;  // WIMP-nucleon cross-sections


// UTILITY FUNCTIONS DECLARATIONS --------------------------------------

/* Provide prototypes for utility functions used by this example
   program.  The function definitions are given later in this file,
   after the main() routine. */

// Print description of input
void WriteDescription(const int type);

// Get WIMP mass and couplings
bool GetWIMPParams(const int type, double& M, double& xpSI, double& xnSI,
                   double& xpSD, double& xnSD);


// MAIN PROGRAM --------------------------------------------------------

int main(int argc, char* argv[]){
  
  int type;
  double M,xpSI,xnSI,xpSD,xnSD,GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an,
         sigmapSI,sigmanSI,sigmapSD,sigmanSD;
  
  // Parse command line options
  // Notably, determining how WIMP parameters will be specified.
  // Default command line option (no argument) will give type = TYPE_MSIGMA.
  type = TYPE_MSIGMA;
  for (int i=1; i<argc; i++) {
    if (string(argv[i]) == "--mG")
      type = TYPE_MG;
    else if (string(argv[i]) == "--mfa")
      type = TYPE_MFA;
    else if (string(argv[i]) == "--msigma")
      type = TYPE_MSIGMA;
    else if (string(argv[i]) == "--help") {
      cout << "Usage:" << endl;
      cout << "  ./DDCalc0_exampleC [--mG|--mfa|--msigma]" << endl;
      cout << "where the optional flag specifies the form in which the WIMP-" << endl;
      cout << "nucleon couplings will be provided (default: --msigma)." << endl;
      exit(0);
    } else {
      cout << "WARNING:  Ignoring unknown argument '" << argv[i] << "'." << endl;
    }
  }
  
  /* Write out directions for specifying input to this example
     program.  WriteDescription is defined above. */
  WriteDescription(type);
  
  /* Initialize the DDCalc0 module. */
  DDCalc0_Init();

  /* Initialize any experiments to be used.  The argument indicates if
     extra sub-interval calculations should be performed.  Those
     calculations are required for maximum gap analyses, but are
     unnecessary for calculating total rates and likelihoods.  If
     'false' is given as an argument, a no-background-subtraction
     p-value can still be calculated, but a Poisson is used instead of
     the maximum gap. We show some maximum gap results below, so we
     must use 'true' here.  The default value for this optional
     argument is 'true', so the explicit argument below is not
     actually necessary here. */
  XENON100_2012_Init(true);
  LUX_2013_Init(true);
  SuperCDMS_2014_Init(true);
  DARWIN_Ar_2015_Init(true);
  DARWIN_Xe_2015_Init(true);

  /* Can optionally specify a minimum recoil energy to be included in
     the rate calculations [keV].  Note the efficiency curves already
     account for detector and analysis thresholds regardless of this
     setting, so setting this to 0 keV (the default behavior when
     initialization is performed) does not imply that very low energy
     recoils actually contribute to the signal. */
  // EXAMPLE: Uncomment to set a minimum recoil energy of 3 keV.
  //XENON100_2012_SetEmin(3.0);
  //LUX_2013_SetEmin(3.0);
  //SuperCDMS_2014_SetEmin(3.0);
  //DARWIN_Ar_2015_SetEmin(3.0);
  //DARWIN_Xe_2015_SetEmin(3.0);

  /* Optionally set the Standard Halo Model parameters:
       rho     Local dark matter density [GeV/cm^3]
       vrot    Local disk rotation speed [km/s]
       v0      Maxwell-Boltzmann most probable speed [km/s]
       vesc    Galactic escape speed [km/s]
     This example uses the default values (and is thus optional). */
  //DDCalc0_SetSHM(0.4, 235.0, 235.0, 550.0)
  
  // INPUT LOOP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Loop over input to this example program.
  // GetWIMPParams is defined below.
  while(GetWIMPParams(type,M,xpSI,xnSI,xpSD,xnSD)) {

    cout << endl;
    
    /* Set the WIMP parameters.
       There are three ways to specify the WIMP-nucleon couplings, with
       the WIMP mass [GeV] always the first argument:
         * DDCalc0_SetWIMP_mfa(m,fp,fn,ap,an)
           The standard couplings fp,fn [GeV^-2] & ap,an [unitless]
         * DDCalc0_SetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD)
           The effective 4 fermion vertex couplings GpSI,GnSI,GpSD,GnSD
           [GeV^-2], related by:
               GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
               GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
         * DDCalc0_SetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
           The WIMP-nucleon cross-sections [pb] (use a negative value
           to indicate the corresponding coupling should be negative).
       In the above, 'p' is for proton, 'n' is for neutron, 'SI' is for
       spin-independent, and 'SD' is for spin-dependent. */
    switch (type) {
      case TYPE_MG:
        DDCalc0_SetWIMP_mG(M,xpSI,xnSI,xpSD,xnSD);
        break;
      case TYPE_MFA:
        DDCalc0_SetWIMP_mfa(M,xpSI,xnSI,xpSD,xnSD);
        break;
      case TYPE_MSIGMA:
        DDCalc0_SetWIMP_msigma(M,xpSI,xnSI,xpSD,xnSD);
        break;
    }
    
    /* Get the WIMP parameters with the same signatures and units as
       above.  The only difference is that WIMP-nucleon cross-sections
       are always positive. */
    DDCalc0_GetWIMP_mfa(M,fp,fn,ap,an);
    DDCalc0_GetWIMP_mG(M,GpSI,GnSI,GpSD,GnSD);
    DDCalc0_GetWIMP_msigma(M,sigmapSI,sigmanSI,sigmapSD,sigmanSD);
    
    /* Print out the above WIMP mass, couplings, and cross sections. */
    printf("%s %- #12.5g\n","WIMP mass [GeV]     ",M);
    cout << endl;
    printf("%-28s %11s %11s %11s %11s\n","WIMP-nucleon couplings",
           " proton-SI "," neutron-SI"," proton-SD "," neutron-SD");
    printf("%-28s %- #11.5g %- #11.5g %- #11.5g %- #11.5g\n",
           "  G [GeV^-2]",GpSI,GnSI,GpSD,GnSD);
    printf("%-28s %- #11.5g %- #11.5g %- #11.5g %- #11.5g\n",
           "  f & a [GeV^-2,unitless]",fp,fn,ap,an);
    printf("%-28s %- #11.5g %- #11.5g %- #11.5g %- #11.5g\n",
           "  cross-section [pb]",sigmapSI,sigmanSI,sigmapSD,sigmanSD);
    cout << endl;
    
    /* Do rate calculations.  After any change to the WIMP or halo
       parameters, perform the rate calculations necessary for predicted
       signals, likelihoods, and/or maximum gap statistics. */
    XENON100_2012_CalcRates();
    LUX_2013_CalcRates();
    SuperCDMS_2014_CalcRates();
    DARWIN_Ar_2015_CalcRates();
    DARWIN_Xe_2015_CalcRates();
    
    /* Header */
    printf("%-20s  %11s  %11s  %11s  %11s  %11s\n","",
           " XENON 2012"," LUX 2013  ","SuCDMS 2014"," DARWIN Ar "," DARWIN Xe ");
    //printf("%-20s  %11s  %11s  %11s  %11s  %11s\n","",
    //       "-----------","-----------","-----------","-----------","-----------");
    
    /* Event quantities. */
    printf("%-20s  % 6i       % 6i       % 6i       % 6i       % 6i     \n",
           "Observed events     ",
           XENON100_2012_Events(),LUX_2013_Events(),
           SuperCDMS_2014_Events(),
           DARWIN_Ar_2015_Events(),DARWIN_Xe_2015_Events());
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Expected background ",
           XENON100_2012_Background(),LUX_2013_Background(),
           SuperCDMS_2014_Background(),
           DARWIN_Ar_2015_Background(),DARWIN_Xe_2015_Background());
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Expected signal     ",
           XENON100_2012_Signal(),LUX_2013_Signal(),
           SuperCDMS_2014_Signal(),
           DARWIN_Ar_2015_Signal(),DARWIN_Xe_2015_Signal());
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "  spin-independent  ",
           XENON100_2012_SignalSI(),LUX_2013_SignalSI(),
           SuperCDMS_2014_SignalSI(),
           DARWIN_Ar_2015_SignalSI(),DARWIN_Xe_2015_SignalSI());
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "  spin-dependent    ",
           XENON100_2012_SignalSD(),LUX_2013_SignalSD(),
           SuperCDMS_2014_SignalSD(),
           DARWIN_Ar_2015_SignalSD(),DARWIN_Xe_2015_SignalSD());
    
    /* The log-likelihoods for the current WIMP; note these are _not_
       multiplied by -2.  The likelihood is calculated using a Poisson
       given the observed x.length() of events and expected signal +
       background. */
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Log-likelihood      ",
           XENON100_2012_LogLikelihood(),LUX_2013_LogLikelihood(),
           SuperCDMS_2014_LogLikelihood(),
           DARWIN_Ar_2015_LogLikelihood(),DARWIN_Xe_2015_LogLikelihood());
    
    /* The logarithm of the p-value, calculated without background
       subtraction, using either the maximum gap statistic or a Poisson
       statistic, depending on how the detector was initialized.  Note
       that this is actually a conservative upper _bound_ on the p-value
       in the event of an unknown background and is useful for excluding
       WIMP parameters.  However, since it is not a true p-value, it
       should not be interpreted as being related to any particular
       likelihood. */
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Max gap log(p-value)",
           XENON100_2012_LogPValue(),LUX_2013_LogPValue(),
           SuperCDMS_2014_LogPValue(),
           DARWIN_Ar_2015_LogPValue(),DARWIN_Xe_2015_LogPValue());
    
    /* Returns a factor x by which the current WIMP cross-sections must
       be multiplied (sigma -> x*sigma, applied to all four WIMP-nucleon
       cross-sections) to achieve the given p-value (specified by its
       logarithm).  Useful for finding the no-background-subtraction
       exclusion limits.  For example, if setWIMP_msigma(100.0,10.0,
       10.0,0.0,0.0) is called, then x*(10. pb) would be the SI
       cross-section at a WIMP mass of 100 GeV at which the experiment
       is excluded at the 90% CL (p=1-CL). */
    double lnp = log(0.1);  // default value for optional argument
    printf("%-20s  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g  %- #11.5g\n",
           "Max gap x for 90% CL",
           XENON100_2012_ScaleToPValue(lnp),LUX_2013_ScaleToPValue(lnp),
           SuperCDMS_2014_ScaleToPValue(lnp),
           DARWIN_Ar_2015_ScaleToPValue(lnp),DARWIN_Xe_2015_ScaleToPValue(lnp));
    cout << "         * Factor x such that sigma->x*sigma gives desired p-value" << endl;
    //cout << endl;

  }  // END INPUT LOOP <<<<<<<<<<<<<<<<<<<<<<<<<

} 


// UTILITY FUNCTION DEFINITIONS ----------------------------------------

/* Write a description of how input parameters should be specified. */
void WriteDescription(const int type){
  cout << endl;
  cout << "Enter WIMP parameters below.  Only the first two are necessary." << endl;
  cout << "A blank line terminates input.  The parameters are:" << endl;
  //cout << "Enter WIMP parameters below.  The parameters are:" << endl;
  cout << endl;
  switch (type) {
    case TYPE_MG:
      cout << "  M     WIMP mass [GeV]" << endl;
      cout << "  GpSI  Spin-independent WIMP-proton effective coupling [GeV^-2]" << endl;
      cout << "  GnSI  Spin-independent WIMP-neutron effective coupling [GeV^-2]" << endl;
      cout << "  GpSD  Spin-dependent WIMP-proton effective coupling [GeV^-2]" << endl;
      cout << "  GnSD  Spin-dependent WIMP-neutron effective coupling [GeV^-2]" << endl;
      break;
    case TYPE_MFA:
      cout << "  M     WIMP mass [GeV]" << endl;
      cout << "  fp    Spin-independent WIMP-proton effective coupling [GeV^-2]" << endl;
      cout << "  fn    Spin-independent WIMP-neutron effective coupling [GeV^-2]" << endl;
      cout << "  ap    Spin-dependent WIMP-proton effective coupling [unitless]" << endl;
      cout << "  an    Spin-dependent WIMP-neutron effective coupling [unitless]" << endl;
      break;
    case TYPE_MSIGMA:
      cout << "  M         WIMP mass [GeV]" << endl;
      cout << "  sigmapSI  Spin-independent WIMP-proton cross-section [pb]" << endl;
      cout << "  sigmanSI  Spin-independent WIMP-neutron cross-section [pb]" << endl;
      cout << "  sigmapSD  Spin-dependent WIMP-proton cross-section [pb]" << endl;
      cout << "  sigmanSD  Spin-dependent WIMP-neutron cross-section [pb]" << endl;
      cout << endl;
      cout << "Negative cross-section values can be given to indicate the" << endl;
      cout << "corresponding coupling should be taken to be negative." << endl;
      break;
  }
}


/* Read WIMP parameters (mass & couplings) from standard input. */
bool GetWIMPParams(const int type, double& M, double& xpSI, double& xnSI,
                   double& xpSD, double& xnSD) {
  cout << endl;
  cout << "------------------------------------------------------------" << endl;
  switch (type) {
    case TYPE_MG:
      cout << "Enter values <M GpSI GnSI GpSD GnSD>:" << endl;
      break;
    case TYPE_MFA:
      cout << "Enter values <M fp fn ap an>:" << endl;
      break;
    case TYPE_MSIGMA:
      cout << "Enter values <M sigmapSI sigmanSI sigmapSD sigmanSD>:" << endl;
      break;
    }
  
  // Get input line for parsing
  string line;
  getline(cin, line);
  istringstream iss(line);

  // Parse input line
  if (!(iss >> M)) return false;
  if (!(iss >> xpSI)) return false;
  if (!(iss >> xnSI)) {
    xnSI = xpSI; xpSD = 0.0; xnSD = 0.0;
    return true;
  }
  if (!(iss >> xpSD)) {
    xpSD = 0.0; xnSD = 0.0;
    return true;
  }
  if (!(iss >> xnSD)) {
    xnSD = xpSD;
    return true;
  }
  return true;
}


// END FILE ------------------------------------------------------------

 
