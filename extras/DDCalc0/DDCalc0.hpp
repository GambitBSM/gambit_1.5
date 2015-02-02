/**********************************************************************
 *
 * Interface to various DDCalc0 routines.
 *
 *      C Savage       Nordita          2014
 *      A Scaffidi     U of Adelaide    2014
 *      M White        U of Adelaide    2014
 *
 **********************************************************************/

#ifndef DDCalc0_H
#define DDCalc0_H

// INTERNAL ROUTINES ---------------------------------------------------
// Declare routines found in DDCalc0 module library.
// Not intended to be called directly by users.  If you want to know how
// to call DDCalc0, see the EXTERNAL ROUTINES section below.

// Naming conventions for the object/library symbols of Fortran module
// routines are typically:
//   __<modulename>_MOD_<routinename>  [gfortran]
//   <modulename>_mp_<routinename>_    [ifort]
// where the module and routine names are in lower case.  To avoid
// compiler-dependendent symbol names, BIND() statements are used in
// the Fortran source code to explicitly specify the symbol names.
// We chose as our naming convention the following:
//   C_DDCALC0_<routinename>
// where the routine name is in lower case.  We call wrapper versions
// of the routines that are C++ type compatible (hence the 'C_' prefix).
// [n.b.: For non-module routines, '<routinename>_' (again in lower
// case) is the convention for both compilers.]
// 
// Arguments are explained later.
// 
extern "C" {
  // Initialization
  void C_DDCALC0_ddcalc0_init();
  
  // Experiment initialization
  void C_DDCALC0_xenon100_2012_init(const bool *);
  void C_DDCALC0_lux_2013_init(const bool *);
  void C_DDCALC0_darwin_ar_2015_init(const bool *);
  void C_DDCALC0_darwin_xe_2015_init(const bool *);
  
  // Halo
  void C_DDCALC0_ddcalc0_setshm(const double *,
                  const double *, const double *, const double *);
  
  // WIMP parameters (set/get)
  void C_DDCALC0_ddcalc0_setwimp_mfa(const double *,
                  const double *, const double *, const double *, const double *);
  void C_DDCALC0_ddcalc0_setwimp_mg(const double *,
                  const double *, const double *, const double *, const double *);
  void C_DDCALC0_ddcalc0_setwimp_msigma(const double *,
                  const double *, const double *, const double *, const double *);
  void C_DDCALC0_ddcalc0_getwimp_mfa(double *,
                  double *, double *, double *, double *);
  void C_DDCALC0_ddcalc0_getwimp_mg(double *,
                  double *, double *, double *, double *);
  void C_DDCALC0_ddcalc0_getwimp_msigma(double *,
                  double *, double *, double *, double *);
  
  // Minimum recoil energy to consider
  void C_DDCALC0_xenon100_2012_setemin(const double *);
  void C_DDCALC0_lux_2013_setemin(const double *);
  void C_DDCALC0_darwin_ar_2015_setemin(const double *);
  void C_DDCALC0_darwin_xe_2015_setemin(const double *);
  
  // Perform calculations
  void C_DDCALC0_xenon100_2012_calcrates();
  void C_DDCALC0_lux_2013_calcrates();
  void C_DDCALC0_darwin_ar_2015_calcrates();
  void C_DDCALC0_darwin_xe_2015_calcrates();
  
  // Number of events
  int C_DDCALC0_xenon100_2012_events();
  int C_DDCALC0_lux_2013_events();
  int C_DDCALC0_darwin_ar_2015_events();
  int C_DDCALC0_darwin_xe_2015_events();
  
  // Expected backgrounds
  double C_DDCALC0_xenon100_2012_background();
  double C_DDCALC0_lux_2013_background();
  double C_DDCALC0_darwin_ar_2015_background();
  double C_DDCALC0_darwin_xe_2015_background();
  
  // Expected signal
  double C_DDCALC0_xenon100_2012_signal();
  double C_DDCALC0_lux_2013_signal();
  double C_DDCALC0_darwin_ar_2015_signal();
  double C_DDCALC0_darwin_xe_2015_signal();
  
  // Expected signal (spin-independent)
  double C_DDCALC0_xenon100_2012_signalsi();
  double C_DDCALC0_lux_2013_signalsi();
  double C_DDCALC0_darwin_ar_2015_signalsi();
  double C_DDCALC0_darwin_xe_2015_signalsi();
  
  // Expected signal (spin-dependent)
  double C_DDCALC0_xenon100_2012_signalsd();
  double C_DDCALC0_lux_2013_signalsd();
  double C_DDCALC0_darwin_ar_2015_signalsd();
  double C_DDCALC0_darwin_xe_2015_signalsd();
  
  // Log-likelihood
  double C_DDCALC0_xenon100_2012_loglikelihood();
  double C_DDCALC0_lux_2013_loglikelihood();
  double C_DDCALC0_darwin_ar_2015_loglikelihood();
  double C_DDCALC0_darwin_xe_2015_loglikelihood();
  
  // Log of the p-value
  double C_DDCALC0_xenon100_2012_logpvalue();
  double C_DDCALC0_lux_2013_logpvalue();
  double C_DDCALC0_darwin_ar_2015_logpvalue();
  double C_DDCALC0_darwin_xe_2015_logpvalue();
  
  // Factor x by which sigma -> x*sigma would yield
  // given p-value (given as log(p))
  double C_DDCALC0_xenon100_2012_scaletopvalue(const double *);
  double C_DDCALC0_lux_2013_scaletopvalue(const double *);
  double C_DDCALC0_darwin_ar_2015_scaletopvalue(const double *);
  double C_DDCALC0_darwin_xe_2015_scaletopvalue(const double *);
  
}


// EXTERNAL (USER) ROUTINES --------------------------------------------

// Initialize the module
void DDCalc0_Init() {
  C_DDCALC0_ddcalc0_init();
};

// Initialize any experiments for which likelihoods are to be calculated.
// The argument indicates if extra calculations necessary for maximum gap
// statistics should be performed (unnecessary for likelihoods).
void XENON100_2012_Init(const bool intervals=true) {
  C_DDCALC0_xenon100_2012_init(&intervals);
};
void LUX_2013_Init(const bool intervals=true) {
  C_DDCALC0_lux_2013_init(&intervals);
};
void DARWIN_Ar_2015_Init(const bool intervals=true) {
  C_DDCALC0_darwin_ar_2015_init(&intervals);
};
void DARWIN_Xe_2015_Init(const bool intervals=true) {
  C_DDCALC0_darwin_xe_2015_init(&intervals);
};

// Optionally set the Standard Halo Model parameters:
//    rho     Local dark matter density [GeV/cm^3]
//    vrot    Local disk rotation speed [km/s]
//    v0      Maxwell-Boltzmann most probable speed [km/s]
//    vesc    Galactic escape speed [km/s]
// This example uses the default values (and is thus optional).
void DDCalc0_SetSHM(const double rho=0.4, const double vrot=235.0,
                    const double v0=235.0, const double vesc=550.0) {
  C_DDCALC0_ddcalc0_setshm(&rho,&vrot,&v0,&vesc);
}

// Set the WIMP parameters.
// There are three ways to specify the WIMP-nucleon couplings, with the
// WIMP mass [GeV] always the first argument:
//   * SetWIMP_mfa(m,fp,fn,ap,an)
//     The standard couplings fp,fn [GeV^-2] & ap,an [unitless]
//   * SetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD)
//     The effective 4 fermion vertex couplings GpSI,GnSI,GpSD,GnSD
//     [GeV^-2], related by:
//         GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
//         GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
//   * SetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
//     The WIMP-nucleon cross-sections [pb] (use a negative value
//     to indicate the corresponding coupling should be negative).
// In the above, 'p' is for proton, 'n' is for neutron, 'SI' is for
// spin-independent, and 'SD' is for spin-dependent.
void DDCalc0_SetWIMP_mfa(const double m, const double fp, const double fn,
                         const double ap, const double an) {
  C_DDCALC0_ddcalc0_setwimp_mfa(&m,&fp,&fn,&ap,&an);
}
void DDCalc0_SetWIMP_mG(const double m, const double GpSI, const double GnSI,
                        const double GpSD, const double GnSD) {
  C_DDCALC0_ddcalc0_setwimp_mg(&m,&GpSI,&GnSI,&GpSD,&GnSD);
}
void DDCalc0_SetWIMP_msigma(const double m,
                            const double sigmapSI, const double sigmanSI,
                            const double sigmapSD, const double sigmanSD) {
  C_DDCALC0_ddcalc0_setwimp_msigma(&m,&sigmapSI,&sigmanSI,&sigmapSD,&sigmanSD);
}

// Get the WIMP parameters with the same signatures and units as above.
// The only difference is that WIMP-nucleon cross-sections are always
// positive.
void DDCalc0_GetWIMP_mfa(double& m, double& fp, double& fn,
                         double& ap, double& an) {
  C_DDCALC0_ddcalc0_getwimp_mfa(&m,&fp,&fn,&ap,&an);
}
void DDCalc0_GetWIMP_mG(double& m, double& GpSI, double& GnSI,
                        double& GpSD, double& GnSD) {
  C_DDCALC0_ddcalc0_getwimp_mg(&m,&GpSI,&GnSI,&GpSD,&GnSD);
}
void DDCalc0_GetWIMP_msigma(double& m,
                            double& sigmapSI, double& sigmanSI,
                            double& sigmapSD, double& sigmanSD) {
  C_DDCALC0_ddcalc0_getwimp_msigma(&m,&sigmapSI,&sigmanSI,&sigmapSD,&sigmanSD);
}

// Specify the minimum recoil energy to be included in the rate
// calculations [keV].  Note the efficiency curves already account for
// detector and analysis thresholds regardless of this setting, so
// setting this to 0 keV (the default behavior when initialization is
// performed) does not imply that very low energy recoils actually
// contribute to the signal.
void XENON100_2012_SetEmin(const double Emin) {
  C_DDCALC0_xenon100_2012_setemin(&Emin);
};
void LUX_2013_SetEmin(const double Emin) {
  C_DDCALC0_lux_2013_setemin(&Emin);
};
void DARWIN_Ar_2015_SetEmin(const double Emin) {
  C_DDCALC0_darwin_ar_2015_setemin(&Emin);
};
void DARWIN_Xe_2015_SetEmin(const double Emin) {
  C_DDCALC0_darwin_xe_2015_setemin(&Emin);
};

// After any change to the WIMP or halo parameters, perform the
// rate calculations necessary for the likelihoods.
void XENON100_2012_CalcRates() {
  C_DDCALC0_xenon100_2012_calcrates();
}
void LUX_2013_CalcRates() {
  C_DDCALC0_lux_2013_calcrates();
}
void DARWIN_Ar_2015_CalcRates() {
  C_DDCALC0_darwin_ar_2015_calcrates();
}
void DARWIN_Xe_2015_CalcRates() {
  C_DDCALC0_darwin_xe_2015_calcrates();
}

// Observed events
int XENON100_2012_Events() {
  return C_DDCALC0_xenon100_2012_events();
}
int LUX_2013_Events() {
  return C_DDCALC0_lux_2013_events();
}
int DARWIN_Ar_2015_Events() {
  return C_DDCALC0_darwin_ar_2015_events();
}
int DARWIN_Xe_2015_Events() {
  return C_DDCALC0_darwin_xe_2015_events();
}

// Expected background
double XENON100_2012_Background() {
  return C_DDCALC0_xenon100_2012_background();
}
double LUX_2013_Background() {
  return C_DDCALC0_lux_2013_background();
}
double DARWIN_Ar_2015_Background() {
  return C_DDCALC0_darwin_ar_2015_background();
}
double DARWIN_Xe_2015_Background() {
  return C_DDCALC0_darwin_xe_2015_background();
}

// Expected signal
double XENON100_2012_Signal() {
  return C_DDCALC0_xenon100_2012_signal();
}
double LUX_2013_Signal() {
  return C_DDCALC0_lux_2013_signal();
}
double DARWIN_Ar_2015_Signal() {
  return C_DDCALC0_darwin_ar_2015_signal();
}
double DARWIN_Xe_2015_Signal() {
  return C_DDCALC0_darwin_xe_2015_signal();
}

// Expected signal (spin-independent)
double XENON100_2012_SignalSI() {
  return C_DDCALC0_xenon100_2012_signalsi();
}
double LUX_2013_SignalSI() {
  return C_DDCALC0_lux_2013_signalsi();
}
double DARWIN_Ar_2015_SignalSI() {
  return C_DDCALC0_darwin_ar_2015_signalsi();
}
double DARWIN_Xe_2015_SignalSI() {
  return C_DDCALC0_darwin_xe_2015_signalsi();
}

// Expected signal (spin-dependent)
double XENON100_2012_SignalSD() {
  return C_DDCALC0_xenon100_2012_signalsd();
}
double LUX_2013_SignalSD() {
  return C_DDCALC0_lux_2013_signalsd();
}
double DARWIN_Ar_2015_SignalSD() {
  return C_DDCALC0_darwin_ar_2015_signalsd();
}
double DARWIN_Xe_2015_SignalSD() {
  return C_DDCALC0_darwin_xe_2015_signalsd();
}

// The log-likelihoods for the current WIMP; note these are _not_
// multiplied by -2.  The likelihood is calculated using a Poisson
// given the observed number of events and expected signal + background.
double XENON100_2012_LogLikelihood() {
  return C_DDCALC0_xenon100_2012_loglikelihood();
}
double LUX_2013_LogLikelihood() {
  return C_DDCALC0_lux_2013_loglikelihood();
}
double DARWIN_Ar_2015_LogLikelihood() {
  return C_DDCALC0_darwin_ar_2015_loglikelihood();
}
double DARWIN_Xe_2015_LogLikelihood() {
  return C_DDCALC0_darwin_xe_2015_loglikelihood();
}

// The logarithm of the p-value, calculated without background
// subtraction, using either the maximum gap statistic or a Poisson
// statistic, depending on how the detector was initialized.  Note that
// this is actually a conservative upper _bound_ on the p-value in the
// event of an unknown background and is useful for excluding WIMP
// parameters.  However, since it is not a true p-value, it should not
// be interpreted as being related to any particular likelihood.
double XENON100_2012_LogPValue() {
  return C_DDCALC0_xenon100_2012_logpvalue();
}
double LUX_2013_LogPValue() {
  return C_DDCALC0_lux_2013_logpvalue();
}
double DARWIN_Ar_2015_LogPValue() {
  return C_DDCALC0_darwin_ar_2015_logpvalue();
}
double DARWIN_Xe_2015_LogPValue() {
  return C_DDCALC0_darwin_xe_2015_logpvalue();
}

// Returns a factor x by which the current WIMP cross-sections must
// be multiplied (sigma -> x*sigma, applied to all four WIMP-nucleon
// cross-sections) to achieve the given p-value (specified by its
// logarithm).  Useful for finding the no-background-subtraction
// exclusion limits.  For example, if setWIMP_msigma(100.0,10.0,10.0,
// 0.0,0.0) is called, then x*(10. pb) would be the SI cross-section
// at a WIMP mass of 100 GeV at which the experiment is excluded at
// the 90% CL (p=1-CL).
double XENON100_2012_ScaleToPValue(const double logp=-2.302585) {
  return C_DDCALC0_xenon100_2012_scaletopvalue(&logp);
}
double LUX_2013_ScaleToPValue(const double logp=-2.302585) {
  return C_DDCALC0_lux_2013_scaletopvalue(&logp);
}
double DARWIN_Ar_2015_ScaleToPValue(const double logp=-2.302585) {
  return C_DDCALC0_darwin_ar_2015_scaletopvalue(&logp);
}
double DARWIN_Xe_2015_ScaleToPValue(const double logp=-2.302585) {
  return C_DDCALC0_darwin_xe_2015_scaletopvalue(&logp);
}


#endif  // DDCalc0_H

