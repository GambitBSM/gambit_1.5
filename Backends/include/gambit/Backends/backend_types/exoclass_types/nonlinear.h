#include "common.h"

enum non_linear_method {nl_none,nl_halofit};

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialized by nonlinear_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */

struct nonlinear {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

  enum non_linear_method method; /**< method for computing non-linear corrections (none, Halogit, etc.) */

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

  int pk_size;     /**< k_size = total number of pk: 1 (P_m) if no massive neutrinos, 2 (P_m and P_cb) if massive neutrinos are present*/
  int index_pk_m;
  int index_pk_cb;
  short has_pk_cb; /** calculate P(k) with only cold dark matter and baryons*/
  int k_size;      /**< k_size = total number of k values */
  double * k;      /**< k[index_k] = list of k values */
  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values */

  double ** nl_corr_density;   /**< nl_corr_density[index_pk][index_tau * ppt->k_size + index_k] */
  double ** k_nl;  /**< wavenumber at which non-linear corrections become important, defined differently by different non_linear_method's */
  int index_tau_min_nl; /**< index of smallest value of tau at which nonlinear corrections have been computed (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - parameters for the pk_eq method */

  short has_pk_eq;               /**< flag: will we use the pk_eq method? */

  int index_pk_eq_w;                /**< index of w in table pk_eq_w_and_Omega */
  int index_pk_eq_Omega_m;          /**< index of Omega_m in table pk_eq_w_and_Omega */
  int pk_eq_size;                   /**< number of indices in table pk_eq_w_and_Omega */

  int pk_eq_tau_size;               /**< number of times (and raws in table pk_eq_w_and_Omega) */

  double * pk_eq_tau;               /**< table of time values */
  double * pk_eq_w_and_Omega;       /**< table of background quantites */
  double * pk_eq_ddw_and_ddOmega;   /**< table of second derivatives */

  //@{

  //@}

  /** @name - technical parameters */

  //@{

  short nonlinear_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  //@}
};
