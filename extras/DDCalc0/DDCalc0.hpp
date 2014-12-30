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
// Not intended to be called directly by users (see following section).

// Naming convention for Fortran module routine is typically:
//   __<modulename>_MOD_<routinename>
// where the module and routine names are in lower case.
// Arguments are explained later.
extern "C" {
  // Initialization
  void __ddcalc0_MOD_ddcalc0_init();
  
  // Experiment initialization
  void __ddcalc0_MOD_xenon100_2012_init(const bool *);
  void __ddcalc0_MOD_lux_2013_init(const bool *);
  void __ddcalc0_MOD_darwin_ar_2014_init(const bool *);
  void __ddcalc0_MOD_darwin_xe_2014_init(const bool *);
  
  // Halo
  void __ddcalc0_MOD_ddcalc0_setshm(const double *,
                  const double *, const double *, const double *);
  
  // WIMP parameters (set/get)
  void __ddcalc0_MOD_ddcalc0_setwimp_mfa(const double *,
                  const double *, const double *, const double *, const double *);
  void __ddcalc0_MOD_ddcalc0_setwimp_mg(const double *,
                  const double *, const double *, const double *, const double *);
  void __ddcalc0_MOD_ddcalc0_setwimp_msigma(const double *,
                  const double *, const double *, const double *, const double *);
  void __ddcalc0_MOD_ddcalc0_getwimp_mfa(double *,
                  double *, double *, double *, double *);
  void __ddcalc0_MOD_ddcalc0_getwimp_mg(double *,
                  double *, double *, double *, double *);
  void __ddcalc0_MOD_ddcalc0_getwimp_msigma(double *,
                  double *, double *, double *, double *);
  
  // Minimum recoil energy to consider
  void __ddcalc0_MOD_xenon100_2012_setemin(const double *);
  void __ddcalc0_MOD_lux_2013_setemin(const double *);
  void __ddcalc0_MOD_darwin_ar_2014_setemin(const double *);
  void __ddcalc0_MOD_darwin_xe_2014_setemin(const double *);
  
  // Perform calculations
  void __ddcalc0_MOD_xenon100_2012_calcrates();
  void __ddcalc0_MOD_lux_2013_calcrates();
  void __ddcalc0_MOD_darwin_ar_2014_calcrates();
  void __ddcalc0_MOD_darwin_xe_2014_calcrates();
  
  // Number of events
  int __ddcalc0_MOD_xenon100_2012_events();
  int __ddcalc0_MOD_lux_2013_events();
  int __ddcalc0_MOD_darwin_ar_2014_events();
  int __ddcalc0_MOD_darwin_xe_2014_events();
  
  // Expected backgrounds
  double __ddcalc0_MOD_xenon100_2012_background();
  double __ddcalc0_MOD_lux_2013_background();
  double __ddcalc0_MOD_darwin_ar_2014_background();
  double __ddcalc0_MOD_darwin_xe_2014_background();
  
  // Expected signal
  double __ddcalc0_MOD_xenon100_2012_signal();
  double __ddcalc0_MOD_lux_2013_signal();
  double __ddcalc0_MOD_darwin_ar_2014_signal();
  double __ddcalc0_MOD_darwin_xe_2014_signal();
  
  // Expected signal (spin-independent)
  double __ddcalc0_MOD_xenon100_2012_signalsi();
  double __ddcalc0_MOD_lux_2013_signalsi();
  double __ddcalc0_MOD_darwin_ar_2014_signalsi();
  double __ddcalc0_MOD_darwin_xe_2014_signalsi();
  
  // Expected signal (spin-dependent)
  double __ddcalc0_MOD_xenon100_2012_signalsd();
  double __ddcalc0_MOD_lux_2013_signalsd();
  double __ddcalc0_MOD_darwin_ar_2014_signalsd();
  double __ddcalc0_MOD_darwin_xe_2014_signalsd();
  
  // Log-likelihood
  double __ddcalc0_MOD_xenon100_2012_loglikelihood();
  double __ddcalc0_MOD_lux_2013_loglikelihood();
  double __ddcalc0_MOD_darwin_ar_2014_loglikelihood();
  double __ddcalc0_MOD_darwin_xe_2014_loglikelihood();
  
  // Log of the p-value
  double __ddcalc0_MOD_xenon100_2012_logpvalue();
  double __ddcalc0_MOD_lux_2013_logpvalue();
  double __ddcalc0_MOD_darwin_ar_2014_logpvalue();
  double __ddcalc0_MOD_darwin_xe_2014_logpvalue();
  
  // Factor x by which sigma -> x*sigma would yield
  // given p-value (given as log(p))
  double __ddcalc0_MOD_xenon100_2012_scaletopvalue(double *);
  double __ddcalc0_MOD_lux_2013_scaletopvalue(double *);
  double __ddcalc0_MOD_darwin_ar_2014_scaletopvalue(double *);
  double __ddcalc0_MOD_darwin_xe_2014_scaletopvalue(double *);
  
}


// EXTERNAL (USER) ROUTINES --------------------------------------------

// Initialize the module
void DDCalc0_Init() {
  __ddcalc0_MOD_ddcalc0_init();
};

// Initialize any experiments for which likelihoods are to be calculated.
// The argument indicates if extra calculations necessary for maximum gap
// statistics should be performed (unnecessary for likelihoods).
void XENON100_2012_Init(const bool intervals=true) {
  __ddcalc0_MOD_xenon100_2012_init(&intervals);
};
void LUX_2013_Init(const bool intervals=true) {
  __ddcalc0_MOD_lux_2013_init(&intervals);
};
void DARWIN_Ar_2014_Init(const bool intervals=true) {
  __ddcalc0_MOD_darwin_ar_2014_init(&intervals);
};
void DARWIN_Xe_2014_Init(const bool intervals=true) {
  __ddcalc0_MOD_darwin_xe_2014_init(&intervals);
};

// Optionally set the Standard Halo Model parameters:
//    rho     Local dark matter density [GeV/cm^3]
//    vrot    Local disk rotation speed [km/s]
//    v0      Maxwell-Boltzmann most probable speed [km/s]
//    vesc    Galactic escape speed [km/s]
// This example uses the default values (and is thus optional).
void DDCalc0_SetSHM(const double rho=0.4, const double vrot=235.0,
                    const double v0=235.0, const double vesc=550.0) {
  __ddcalc0_MOD_ddcalc0_setshm(&rho,&vrot,&v0,&vesc);
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
  __ddcalc0_MOD_ddcalc0_setwimp_mfa(&m,&fp,&fn,&ap,&an);
}
void DDCalc0_SetWIMP_mG(const double m, const double GpSI, const double GnSI,
                        const double GpSD, const double GnSD) {
  __ddcalc0_MOD_ddcalc0_setwimp_mg(&m,&GpSI,&GnSI,&GpSD,&GnSD);
}
void DDCalc0_SetWIMP_msigma(const double m,
                            const double sigmapSI, const double sigmanSI,
                            const double sigmapSD, const double sigmanSD) {
  __ddcalc0_MOD_ddcalc0_setwimp_msigma(&m,&sigmapSI,&sigmanSI,&sigmapSD,&sigmanSD);
}

// Get the WIMP parameters with the same signatures and units as above.
// The only difference is that WIMP-nucleon cross-sections are always
// positive.
void DDCalc0_GetWIMP_mfa(double m, double fp, double fn,
                         double ap, double an) {
  __ddcalc0_MOD_ddcalc0_getwimp_mfa(&m,&fp,&fn,&ap,&an);
}
void DDCalc0_GetWIMP_mG(double m, double GpSI, double GnSI,
                        double GpSD, double GnSD) {
  __ddcalc0_MOD_ddcalc0_getwimp_mg(&m,&GpSI,&GnSI,&GpSD,&GnSD);
}
void DDCalc0_GetWIMP_msigma(double m,
                            double sigmapSI, double sigmanSI,
                            double sigmapSD, double sigmanSD) {
  __ddcalc0_MOD_ddcalc0_getwimp_msigma(&m,&sigmapSI,&sigmanSI,&sigmapSD,&sigmanSD);
}

// Specify the minimum recoil energy to be included in the rate
// calculations [keV].  Note the efficiency curves already account for
// detector and analysis thresholds regardless of this setting, so
//  setting this to 0 keV (the default behavior when initialization is
// performed) does not imply that very low energy recoils actually
// contribute to the signal.
void XENON100_2012_SetEmin(const double Emin) {
  __ddcalc0_MOD_ddcalc0_xenon100_2012_setemin(&Emin);
};
void LUX_2013_SetEmin(const double Emin) {
  __ddcalc0_MOD_ddcalc0_lux_2013_setemin(&Emin);
};
void DARWIN_Ar_2014_SetEmin(const double Emin) {
  __ddcalc0_MOD_ddcalc0_darwin_ar_2014_setemin(&Emin);
};
void DARWIN_Xe_2014_SetEmin(const double Emin) {
  __ddcalc0_MOD_ddcalc0_darwin_xe_2014_setemin(&Emin);
};

// After any change to the WIMP or halo parameters, perform the
// rate calculations necessary for the likelihoods.
void XENON100_2012_CalcRates() {
  __ddcalc0_MOD_xenon100_2012_calcrates();
}
void LUX_2013_CalcRates() {
  __ddcalc0_MOD_lux_2013_calcrates();
}
void DARWIN_Ar_2014_CalcRates() {
  __ddcalc0_MOD_darwin_ar_2014_calcrates();
}
void DARWIN_Xe_2014_CalcRates() {
  __ddcalc0_MOD_darwin_xe_2014_calcrates();
}

// Observed events
int XENON100_2012_Events() {
  return __ddcalc0_MOD_xenon100_2012_events();
}
int LUX_2013_Events() {
  return __ddcalc0_MOD_lux_2013_events();
}
int DARWIN_Ar_2014_Events() {
  return __ddcalc0_MOD_darwin_ar_2014_events();
}
int DARWIN_Xe_2014_Events() {
  return __ddcalc0_MOD_darwin_xe_2014_events();
}

// Expected background
double XENON100_2012_Background() {
  return __ddcalc0_MOD_xenon100_2012_background();
}
double LUX_2013_Background() {
  return __ddcalc0_MOD_lux_2013_background();
}
double DARWIN_Ar_2014_Background() {
  return __ddcalc0_MOD_darwin_ar_2014_background();
}
double DARWIN_Xe_2014_Background() {
  return __ddcalc0_MOD_darwin_xe_2014_background();
}

// Expected signal
double XENON100_2012_Signal() {
  return __ddcalc0_MOD_xenon100_2012_signal();
}
double LUX_2013_Signal() {
  return __ddcalc0_MOD_lux_2013_signal();
}
double DARWIN_Ar_2014_Signal() {
  return __ddcalc0_MOD_darwin_ar_2014_signal();
}
double DARWIN_Xe_2014_Signal() {
  return __ddcalc0_MOD_darwin_xe_2014_signal();
}

// Expected signal (spin-independent)
double XENON100_2012_SignalSI() {
  return __ddcalc0_MOD_xenon100_2012_signalsi();
}
double LUX_2013_SignalSI() {
  return __ddcalc0_MOD_lux_2013_signalsi();
}
double DARWIN_Ar_2014_SignalSI() {
  return __ddcalc0_MOD_darwin_ar_2014_signalsi();
}
double DARWIN_Xe_2014_SignalSI() {
  return __ddcalc0_MOD_darwin_xe_2014_signalsi();
}

// Expected signal (spin-dependent)
double XENON100_2012_SignalSD() {
  return __ddcalc0_MOD_xenon100_2012_signalsd();
}
double LUX_2013_SignalSD() {
  return __ddcalc0_MOD_lux_2013_signalsd();
}
double DARWIN_Ar_2014_SignalSD() {
  return __ddcalc0_MOD_darwin_ar_2014_signalsd();
}
double DARWIN_Xe_2014_SignalSD() {
  return __ddcalc0_MOD_darwin_xe_2014_signalsd();
}

// The log-likelihoods for the current WIMP; note these are _not_
// multiplied by -2.  The likelihood is calculated using a Poisson
// given the observed number of events and expected signal + background.
double XENON100_2012_LogLikelihood() {
  return __ddcalc0_MOD_xenon100_2012_loglikelihood();
}
double LUX_2013_LogLikelihood() {
  return __ddcalc0_MOD_lux_2013_loglikelihood();
}
double DARWIN_Ar_2014_LogLikelihood() {
  return __ddcalc0_MOD_darwin_ar_2014_loglikelihood();
}
double DARWIN_Xe_2014_LogLikelihood() {
  return __ddcalc0_MOD_darwin_xe_2014_loglikelihood();
}

// The logarithm of the p-value, calculated without background
// subtraction, using either the maximum gap statistic or a Poisson
// statistic, depending on how the detector was initialized.  Note that
// this is actually a conservative upper _bound_ on the p-value in the
// event of an unknown background and is useful for excluding WIMP
// parameters.  However, since it is not a true p-value, it should not
// be interpreted as being related to any particular likelihood.
double XENON100_2012_LogPValue() {
  return __ddcalc0_MOD_xenon100_2012_logpvalue();
}
double LUX_2013_LogPValue() {
  return __ddcalc0_MOD_lux_2013_logpvalue();
}
double DARWIN_Ar_2014_LogPValue() {
  return __ddcalc0_MOD_darwin_ar_2014_logpvalue();
}
double DARWIN_Xe_2014_LogPValue() {
  return __ddcalc0_MOD_darwin_xe_2014_logpvalue();
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
  return __ddcalc0_MOD_xenon100_2012_scaletopvalue(&logp);
}
double LUX_2013_ScaleToPValue(const double logp=-2.302585) {
  return __ddcalc0_MOD_lux_2013_scaletopvalue(&logp);
}
double DARWIN_Ar_2014_ScaleToPValue(const double logp=-2.302585) {
  return __ddcalc0_MOD_darwin_ar_2014_scaletopvalue(&logp);
}
double DARWIN_Xe_2014_ScaleToPValue(const double logp=-2.302585) {
  return __ddcalc0_MOD_darwin_xe_2014_scaletopvalue(&logp);
}





#endif  // DDCalc0_H

