#include "common.h"

/**
 * List of possible recombination algorithms.
 */

enum recombination_algorithm {
  recfast,
  hyrec
};

/**
 * List of possible reionization schemes.
 */

enum reionization_parametrization {
  reio_none, /**< no reionization */
  reio_camb,  /**< reionization parameterized like in CAMB */
  reio_bins_tanh,  /**< binned reionization history with tanh inteprolation between bins */
  reio_half_tanh,  /**< half a tanh, instead of the full tanh */
  reio_many_tanh,  /**< similar to reio_camb but with more than one tanh */
  reio_inter       /**< linear interpolation between specified points */
};

/**
 * Is the input parameter the reionization redshift or optical depth?
 */

enum reionization_z_or_tau {
  reio_z,  /**< input = redshift */
  reio_tau /**< input = tau */
};

/**
 * All thermodynamics parameters and evolution that other modules need to know.
 *
 * Once initialized by thermodynamics_init(), contains all the
 * necessary information on the thermodynamics, and in particular, a
 * table of thermodynamical quantities as a function of the redshift,
 * used for interpolation in other modules.
 */

struct thermo
{
  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these parameters
   *   and the content of the 'precision' and 'background' structures) */

  //@{

  double YHe;  /**< \f$ Y_{He} \f$: primordial helium fraction */

  enum recombination_algorithm recombination; /**< recombination code */

  enum reionization_parametrization reio_parametrization; /**< reionization scheme */

  enum reionization_z_or_tau reio_z_or_tau; /**< is the input parameter the reionization redshift or optical depth? */

  double tau_reio; /**< if above set to tau, input value of reionization optical depth */

  double z_reio;   /**< if above set to z,   input value of reionization redshift */

  short compute_cb2_derivatives; /**< do we want to include in computation derivatives of baryon sound speed? */

  short compute_damping_scale; /**< do we want to compute the simplest analytic approximation to the photon damping (or diffusion) scale? */

  /** parameters for reio_camb */

  double reionization_width; /**< width of H reionization */

  double reionization_exponent; /**< shape of H reionization */

  double helium_fullreio_redshift; /**< redshift for of helium reionization */

  double helium_fullreio_width; /**< width of helium reionization */

  /** parameters for reio_bins_tanh */

  int binned_reio_num; /**< with how many bins do we want to describe reionization? */

  double * binned_reio_z; /**< central z value for each bin */

  double * binned_reio_xe; /**< imposed \f$ X_e(z)\f$ value at center of each bin */

  double binned_reio_step_sharpness; /**< sharpness of tanh() step interpolating between binned values */

    /** parameters for reio_many_tanh */

  int many_tanh_num; /**< with how many jumps do we want to describe reionization? */

  double * many_tanh_z; /**< central z value for each tanh jump */

  double * many_tanh_xe; /**< imposed \f$ X_e(z)\f$ value at the end of each jump (ie at later times)*/

  double many_tanh_width; /**< sharpness of tanh() steps */

    /** parameters for reio_inter */

  int reio_inter_num; /**< with how many jumps do we want to describe reionization? */

  double * reio_inter_z; /**< discrete z values */

  double * reio_inter_xe; /**< discrete \f$ X_e(z)\f$ values */

  /** parameters for energy injection */

  double annihilation; /** parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */

  short has_on_the_spot; /** flag to specify if we want to use the on-the-spot approximation **/

  double decay; /** parameter describing CDM decay (f/tau, see e.g. 1109.6322)*/

  double annihilation_variation; /** if this parameter is non-zero,
				     the function F(z)=(f <sigma*v> /
				     m_cdm)(z) will be a parabola in
				     log-log scale between zmin and
				     zmax, with a curvature given by
				     annihlation_variation (must be
				     negative), and with a maximum in
				     zmax; it will be constant outside
				     this range */

  double annihilation_z; /** if annihilation_variation is non-zero,
			     this is the value of z at which the
			     parameter annihilation is defined, i.e.
			     F(annihilation_z)=annihilation */

  double annihilation_zmax; /** if annihilation_variation is non-zero,
				redshift above which annihilation rate
				is maximal */

  double annihilation_zmin; /** if annihilation_variation is non-zero,
				redshift below which annihilation rate
				is constant */

  double annihilation_f_halo; /** takes the contribution of DM annihilation in halos into account*/
  double annihilation_z_halo; /** characteristic redshift for DM annihilation in halos*/

  //@}

  /** @name - all indices for the vector of thermodynamical (=th) quantities stored in table */

  //@{

  int index_th_xe;            /**< ionization fraction \f$ x_e \f$ */
  int index_th_dkappa;        /**< Thomson scattering rate \f$ d \kappa / d \tau\f$ (units 1/Mpc) */
  int index_th_tau_d;         /**< Baryon drag optical depth */
  int index_th_ddkappa;       /**< scattering rate derivative \f$ d^2 \kappa / d \tau^2 \f$ */
  int index_th_dddkappa;      /**< scattering rate second derivative \f$ d^3 \kappa / d \tau^3 \f$ */
  int index_th_exp_m_kappa;  /**< \f$ exp^{-\kappa} \f$ */
  int index_th_g;             /**< visibility function \f$ g = (d \kappa / d \tau) * exp^{-\kappa} \f$ */
  int index_th_dg;            /**< visibility function derivative \f$ (d g / d \tau) \f$ */
  int index_th_ddg;           /**< visibility function second derivative \f$ (d^2 g / d \tau^2) \f$ */
  int index_th_Tb;            /**< baryon temperature \f$ T_b \f$ */
  int index_th_cb2;           /**< squared baryon sound speed \f$ c_b^2 \f$ */
  int index_th_dcb2;          /**< derivative wrt conformal time of squared baryon sound speed \f$ d [c_b^2] / d \tau \f$ (only computed if some non-minimal tight-coupling schemes is requested) */
  int index_th_ddcb2;         /**< second derivative wrt conformal time of squared baryon sound speed  \f$ d^2 [c_b^2] / d \tau^2 \f$ (only computed if some non0-minimal tight-coupling schemes is requested) */
  int index_th_rate;          /**< maximum variation rate of \f$ exp^{-\kappa}\f$, g and \f$ (d g / d \tau) \f$, used for computing integration step in perturbation module */
  int index_th_r_d;           /**< simple analytic approximation to the photon comoving damping scale */
  int th_size;                /**< size of thermodynamics vector */

  //@}

  /** @name - thermodynamics interpolation tables */

  //@{

  int tt_size; /**< number of lines (redshift steps) in the tables */
  double * z_table; /**< vector z_table[index_z] with values of redshift (vector of size tt_size) */
  double * thermodynamics_table; /**< table thermodynamics_table[index_z*pth->tt_size+pba->index_th] with all other quantities (array of size th_size*tt_size) */

  //@}

  /** @name - table of their second derivatives, used for spline interpolation */

  //@{

  double * d2thermodynamics_dz2_table; /**< table d2thermodynamics_dz2_table[index_z*pth->tt_size+pba->index_th] with values of \f$ d^2 t_i / dz^2 \f$ (array of size th_size*tt_size) */

  //@}


  /** @name - characteristic quantities like redshift, conformal time and sound horizon at recombination */

  //@{

  double z_rec;   /**< z at which the visibility reaches its maximum (= recombination redshift) */
  double tau_rec; /**< conformal time at which the visibility reaches its maximum (= recombination time) */
  double rs_rec;  /**< comoving sound horizon at recombination */
  double ds_rec;  /**< physical sound horizon at recombination */
  double ra_rec;  /**< conformal angular diameter distance to recombination */
  double da_rec;  /**< physical angular diameter distance to recombination */
  double rd_rec;  /**< comoving photon damping scale at recombination */
  double z_d;     /**< baryon drag redshift */
  double tau_d;   /**< baryon drag time */
  double ds_d;    /**< physical sound horizon at baryon drag */
  double rs_d;    /**< comoving sound horizon at baryon drag */
  double tau_cut; /**< at at which the visibility goes below a fixed fraction of the maximum visibility, used for an approximation in perturbation module */
  double angular_rescaling; /**< [ratio ra_rec / (tau0-tau_rec)]: gives CMB rescaling in angular space relative to flat model (=1 for curvature K=0) */
  double tau_free_streaming;   /**< minimum value of tau at which sfree-streaming approximation can be switched on */

  //@}

  /** @name - initial conformal time at which thermodynamical variables have been be integrated */

  //@{

  double tau_ini; /**< initial conformal time at which thermodynamical variables have been be integrated */

  //@}

/** @name - total number density of electrons today (free or not) */

  //@{

  double n_e; /**< total number density of electrons today (free or not) */

  //@}

  /**
   *@name - some flags needed for thermodynamics functions
   */

  //@{

  short inter_normal;  /**< flag for calling thermodynamics_at_z and find position in interpolation table normally */
  short inter_closeby; /**< flag for calling thermodynamics_at_z and find position in interpolation table starting from previous position in previous call */

  //@}

  /** @name - technical parameters */

  //@{

  short thermodynamics_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}

};
