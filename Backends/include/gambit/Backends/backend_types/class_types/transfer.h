#include "common.h"

/**
 * Structure containing everything about transfer functions in
 * harmonic space \f$ \Delta_l^{X} (q) \f$ that other modules need to
 * know.
 *
 * Once initialized by transfer_init(), contains all tables of
 * transfer functions used for interpolation in other modules, for all
 * requested modes (scalar/vector/tensor), initial conditions, types
 * (temperature, polarization, etc), multipoles l, and wavenumbers q.
 *
 * Wavenumbers are called q in this module and k in the perturbation
 * module. In flat universes k=q. In non-flat universes q and k differ
 * through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector,
 * tensor. q should be used throughout the transfer module, except
 * when interpolating or manipulating the source functions S(k,tau)
 * calculated in the perturbation module: for a given value of q, this
 * should be done at the corresponding k(q).
 *
 * The content of this structure is entirely computed in this module,
 * given the content of the 'precision', 'bessels', 'background',
 * 'thermodynamics' and 'perturbation' structures.
 */

struct transfers {

  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of previous structures) */

  //@{

  double lcmb_rescale; /**< normally set to one, can be used
                          exceptionally to rescale by hand the CMB
                          lensing potential */
  double lcmb_tilt;    /**< normally set to zero, can be used
                          exceptionally to tilt by hand the CMB
                          lensing potential */
  double lcmb_pivot;   /**< if lcmb_tilt non-zero, corresponding pivot
                          scale */

  double selection_bias[_SELECTION_NUM_MAX_];               /**< light-to-mass bias in the transfer function of density number count */
  double selection_magnification_bias[_SELECTION_NUM_MAX_]; /**< magnification bias in the transfer function of density number count */

  short has_nz_file;     /**< Has dN/dz (selection function) input file? */
  short has_nz_analytic; /**< Use analytic form for dN/dz (selection function) distribution? */
  FileName nz_file_name; /**< dN/dz (selection function) input file name */
  int nz_size;           /**< number of redshift values in input tabulated selection function */
  double * nz_z;         /**< redshift values in input tabulated selection function */
  double * nz_nz;        /**< input tabulated values of selection function */
  double * nz_ddnz;      /**< second derivatives in splined selection function*/

  short has_nz_evo_file;      /**< Has dN/dz (evolution function) input file? */
  short has_nz_evo_analytic;  /**< Use analytic form for dN/dz (evolution function) distribution? */
  FileName nz_evo_file_name;  /**< dN/dz (evolution function) input file name */
  int nz_evo_size;            /**< number of redshift values in input tabulated evolution function */
  double * nz_evo_z;          /**< redshift values in input tabulated evolution function */
  double * nz_evo_nz;         /**< input tabulated values of evolution function */
  double * nz_evo_dlog_nz;    /**< log of tabulated values of evolution function */
  double * nz_evo_dd_dlog_nz; /**< second derivatives in splined log of evolution function */

  //@}

  /** @name - flag stating whether we need transfer functions at all */

  //@{

  short has_cls; /**< copy of same flag in perturbation structure */

  //@}

  /** @name - number of modes and transfer function types */

  //@{

  int md_size;       /**< number of modes included in computation */

  int index_tt_t0;      /**< index for transfer type = temperature (j=0 term) */
  int index_tt_t1;      /**< index for transfer type = temperature (j=1 term) */
  int index_tt_t2;      /**< index for transfer type = temperature (j=2 term) */
  int index_tt_e;       /**< index for transfer type = E-polarization */
  int index_tt_b;       /**< index for transfer type = B-polarization */
  int index_tt_lcmb;    /**< index for transfer type = CMB lensing */
  int index_tt_density; /**< index for first bin of transfer type = matter density */
  int index_tt_lensing; /**< index for first bin of transfer type = galaxy lensing */

  int index_tt_rsd;     /**< index for first bin of transfer type = redshift space distortion of number count */
  int index_tt_d0;      /**< index for first bin of transfer type = doppler effect for of number count (j=0 term) */
  int index_tt_d1;      /**< index for first bin of transfer type = doppler effect for of number count (j=1 term) */
  int index_tt_nc_lens; /**< index for first bin of transfer type = lensing for of number count */
  int index_tt_nc_g1;   /**< index for first bin of transfer type = gravity term G1 for of number count */
  int index_tt_nc_g2;   /**< index for first bin of transfer type = gravity term G2 for of number count */
  int index_tt_nc_g3;   /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g4;   /**< index for first bin of transfer type = gravity term G3 for of number count */
  int index_tt_nc_g5;   /**< index for first bin of transfer type = gravity term G3 for of number count */

  int * tt_size;     /**< number of requested transfer types tt_size[index_md] for each mode */

  //@}

  /** @name - number and list of multipoles */

  //@{

  int ** l_size_tt;  /**< number of multipole values for which we effectively compute the transfer function,l_size_tt[index_md][index_tt] */

  int * l_size;   /**< number of multipole values for each requested mode, l_size[index_md] */

  int l_size_max; /**< greatest of all l_size[index_md] */

  int * l;        /**< list of multipole values l[index_l] */

  //int * l_size_bessel; /**< for each wavenumber, maximum value of l at which bessel functions must be evaluated */

  double angular_rescaling; /**< correction between l and k space due to curvature (= comoving angular diameter distance to recombination / comoving radius to recombination) */

  //@}

  /** @name - number and list of wavenumbers */

  //@{

  size_t q_size; /**< number of wavenumber values */

  double * q;  /**< list of wavenumber values, q[index_q] */

  double ** k; /**< list of wavenumber values for each requested mode, k[index_md][index_q]. In flat universes k=q. In non-flat universes q and k differ through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector, tensor. q should be used throughout the transfer module, excepted when interpolating or manipulating the source functions S(k,tau): for a given value of q this should be done in k(q). */

  int index_q_flat_approximation; /**< index of the first q value using the flat rescaling approximation */

  //@}

  /** @name - transfer functions */

  //@{

  double ** transfer; /**< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * ptr->q_size + index_q] */

  //@}

  /** @name - technical parameters */

  //@{

  short initialise_HIS_cache; /**< only true if we are using CLASS for setting up a cache of HIS structures */

  short transfer_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};
