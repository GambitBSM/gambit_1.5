//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  CosmoBit routines relating to axions and 
///  axion-like particles.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@pimperial.ac.uk)
///  \date 2017 Jul
///  \date 2018 May
///  \date 2018 Aug - Sep
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2018 Jan - May
///  \date 2019 Jan, Feb, June, Nov
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 June
///  \date 2019 Mar,June
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2019 June, Nov
///
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///  \date 2020 Mar
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2018 Mar
///  \date 2019 Jul
///  \date 2020 Apr
///
///  *********************************************

#include <valarray>
#include <stdint.h> // save memory addresses as int

#include <gsl/gsl_odeiv2.h>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/CosmoBit/CosmoBit_rollcall.hpp"
#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"

namespace Gambit
{

  namespace CosmoBit
  {
    using namespace LogTags;

    /// Lifetime (in s) of an ALP if only the decay a -> g g is open.
    void lifetime_ALP_agg(double& result)
    {
      using namespace Pipes::lifetime_ALP_agg;

      double gagg = *Param["gagg"]; // in GeV^-1
      double ma = *Param["ma0"]; // in eV

      // Calculate the decay width (in GeV)
      // (It's maybe worth being a separate capability)
      double Gamma = 1/64./pi * pow(gagg,2) * pow(ma,3) * 1e-27;

      // For the lifetime take 1/Gamma and translate GeV^-1 into s by multiplication with "hbar"
      result = 1./Gamma * hbar;

      // Reject points which have a lifetime bigger than 1e17s (or whatever the user chooses)
      // Gets only triggered if the user wishes to do so.
      // !! This is not a real physical bound but it is more to deal with the lack of likelihoods so far. !!
      static bool do_rejection = runOptions->getValueOrDef<bool>(false,"do_rejection");
      static double tdec_max = runOptions->getValueOrDef<double>(1.0e17,"reject_tau_bigger_than");
      if (do_rejection && result > tdec_max)
      {
        std::ostringstream err;
        err << "ALP lifetime (" << result << " [s]) exceeds the threshold of " << tdec_max <<" [s].";
        invalid_point().raise(err.str());
      }
    }

    /// Compute the abundance of ALPs expected from thermal production via Primakoff processes
    void minimum_abundance_ALP(double& result)
    {
      using namespace Pipes::minimum_abundance_ALP;

      double gagg = *Param["gagg"];                 // Read axion-photon coupling in GeV^-1
      double T_R_in_GeV = 1e-3 * (*Param["T_R"]);   // Read reheating temperature in MeV and convert it to GeV

      // Check for stupid input (T_R < m_e) and throw an error if the user really pushed it that far.
      if (m_electron >= T_R_in_GeV)
        CosmoBit_error().raise(LOCAL_INFO,"The reheating temperature is below the electron mass.");

      result = 1.56e-5 * pow(gagg,2) * m_planck * (T_R_in_GeV - m_electron);

    }

    /// Compute the minimal fraction of dark matter in ALPs expected from thermal production via Primakoff processes
    void minimum_fraction_ALP(double& result)
    {
      using namespace Pipes::minimum_fraction_ALP;

      const double rho0_crit_by_h2 = 3.*pow(m_planck_red*1e9,2) * pow((1e5*1e9*hbar/_Mpc_SI_),2); // rho0_crit/(h^2)

      double ma0 = *Param["ma0"];                    // non-thermal ALP mass in eV
      double T = *Param["T_cmb"];                        // CMB temperature in K
      double omega_cdm = *Param["omega_cdm"];        // omega_cdm = Omega_cdm * h^2

      // Consistency check: if the ALP abundance from all thermal processes is less than that expected just from Primakoff processes, invalidate this point.
      double Ya0_min = *Dep::minimum_abundance;
      double ssm0 = CosmoBit_utils::entropy_density_SM(T);           // SM entropy density today in eV^3 (cf. footnote 24 of PDG2018-Astrophysical parameters)
      double rho0_cdm = omega_cdm * rho0_crit_by_h2; // rho0_cdm = Omega_cdm * rho0_crit;
      double rho0_min = Ya0_min * ma0 * ssm0;        // energy density of axions today in eV^4
      result = rho0_min / rho0_cdm;
    }

    /// The fraction of dark matter in decaying ALPs at the time of production
    void DM_fraction_ALP(double& result)
    {
      using namespace Pipes::DM_fraction_ALP;

      const double t_rec = 1e12;                     // Age of the Universe at recombination in s

      double f0_thermal = *Param["f0_thermal"];      // Fraction of DM in ALPs at production, due to thermal production
      double omega_ma = *Dep::RD_oh2;                // Cosmological density of ALPs from vacuum misalignment

      double omega_cdm = *Param["omega_cdm"];        // omega_cdm = Omega_cdm * h^2
      double tau_a = *Dep::lifetime;                 // ALP lifetime in s

      // Consistency check: if the ALP abundance from all thermal processes is less than that expected just from Primakoff processes, invalidate this point.
      double f0_min = *Dep::minimum_fraction;
      if (f0_thermal < f0_min)
      {
        std::ostringstream err;
        err << "The choice of f0_thermal (" << f0_thermal;
        err << ") is in contradiction with the minimum dark matter fraction f0_min (";
        err << f0_min << ") produced via Primakoff processes.";
        invalid_point().raise(err.str());
      }

      // Compute the total fraction of DM in ALPs at production
      double xi_ini = f0_thermal + omega_ma/omega_cdm;

      // Consistency check: invalidate if there are more ALPs than dark matter at the time of recombination (t ~ 1e12s)
      double xi_at_rec = xi_ini * exp(-t_rec/tau_a );
      if (xi_at_rec  > 1.)
      {
        std::ostringstream err;
        err << "ALPs are over-abundant (n_a > n_cdm) at t = 10^12 s. (n_a/n_cdm = "<< xi_at_rec <<")";
        invalid_point().raise(err.str());
      }

      result = xi_ini;
    }

    /// Return the total abundance (Y = n/s) of ALPs
    /// We assume non relativistic ALPs such that rho = n * m
    void total_DM_abundance_ALP(double& result)
    {
      using namespace Pipes::total_DM_abundance_ALP;

      const double rho0_crit_by_h2 = 3.*pow(m_planck_red*1e9,2) * pow((1e5*1e9*hbar/_Mpc_SI_),2); // rho0_crit/(h^2)
      double omega_cdm = *Param["omega_cdm"];  // omega_cdm = Omega_cdm * h^2
      double fraction = *Dep::DM_fraction;
      double rho0_ALP = omega_cdm * fraction* rho0_crit_by_h2; // rho0_cdm = Omega_cdm * rho0_crit

      double ma0 = *Param["ma0"];
      double TCMB = *Param["T_cmb"];
      double ssm0 = CosmoBit_utils::entropy_density_SM(TCMB);

      result = rho0_ALP / (ma0 * ssm0);
    }

    /// Helper function to solve the RHS of the differential equation:
    ///  dT/dt = 15/pi^2 (m_a n_a(t)/ tau_a) T^(-3) - H(T) T , where H(T) = 3.7978719e-7*T*T  
    // @TODO: refer to Eq. number in paper when ready
    ///  params: stores (m_a n_a(t)/ tau_a)
    ///  y[0]: stores SM T[t0]
    int diff_eq_rhs (double t, const double y[], double f[], void *params)
    {
      CosmoBit_utils::fast_interpolation injection_inter = *(static_cast<CosmoBit_utils::fast_interpolation*>(params));
      f[0] = (15.0/(4.0*pi*pi)) * injection_inter.interp(t)/pow(y[0], 3) - 3.7978719e-7*y[0]*y[0]*y[0];
      return GSL_SUCCESS;
    }

    /// @TODO function definition needed
    void compute_dNeff_etaBBN_ALP(map_str_dbl &result)  // takes about 0.2 s with 3000 grid points and 1e-6 convergence criterion atm
    {
      using namespace Pipes::compute_dNeff_etaBBN_ALP;

      double dNeff, Neff_SM, Neff_ALP, eta_ratio;
      double temp_eta_ratio = 1; // temporary values to check convergence of iteration
      double temp_dNeff = 1;
      int ii=0;

      // --- precision parameters --
      double hstart = runOptions->getValueOrDef<double>(1e-6,"hstart"); // initial step size for GSL ODE solver
      double epsrel = runOptions->getValueOrDef<double>(1e-6,"epsrel"); // desired relative accuracy for GSL ODE solver and dNeff & etaBBN
      double max_iter = runOptions->getValueOrDef<int>(10,"max_iter"); // maximum number of iterations before error message is thrown if result not converged
      double grid_size = runOptions->getValueOrDef<int>(3000,"grid_size"); // number of (logarithmic) grid points in t

      double t0 = runOptions->getValueOrDef<double>(1e4,"t_initial"); // initial time in seconds
      double tf = runOptions->getValueOrDef<double>(1e12,"t_final"); // final time in seconds

      std::valarray<double> t_grid(grid_size), T_evo(grid_size), Tnu_grid(grid_size), na_grid(grid_size), injection_grid(grid_size); // arrays containing time dependence of variables

      SM_time_evo SM(t0,tf,grid_size);  // set time evolution of SM
      t_grid = SM.get_t_grid();         // has to be updated when solving the differential equation
      T_evo = SM.get_T_evo();

      // --- model parameters ----
      double Ya0 = *Dep::total_DM_abundance;
      double T0 = T_evo[0];
      double ssm_at_T0 = CosmoBit_utils::entropy_density_SM(T0, true);; // T0 in units of keV, set T_in_eV=True to interpret it correctly

      double na_t0 = Ya0 * ssm_at_T0;     // initial number density of a at t=t0, in units keV^3.
      double m_a = 1e-3*(*Param["ma0"]);  // mass of a in keV
      double tau_a = *Dep::lifetime;      // lifetime of a in seconds

      // loop over iterations until convergence or max_iter is reached
      while( ((fabs(temp_eta_ratio-eta_ratio) > epsrel) || fabs(temp_dNeff - dNeff) > epsrel) && (ii <= max_iter) )
      {
        temp_eta_ratio = eta_ratio;  // to check for convergence
        temp_dNeff = dNeff;

        SM.calc_H_int();
        na_grid = na_t0*exp(-3.0*SM.get_H_int())*exp(-t_grid/tau_a);    // na(t) in units keV^3
        injection_grid = (m_a/tau_a)*na_grid;                           // m_a*na(t)/tau_a in units keV^4/s
        Tnu_grid = SM.get_Tnu_evo()[0]*exp(-1.0*SM.get_H_int());        // T_nu(t) in units keV
        CosmoBit_utils::fast_interpolation injection_inter(t_grid, injection_grid);     // interpolating function

        gsl_odeiv2_system sys = {diff_eq_rhs, NULL, 1, &injection_inter};
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, hstart, 0.0, epsrel);

        double y[1] = {T_evo[0]}; // initial condition: T(t0) = T_SM(t0)
        double t=t0;

        for (int jj = 0; jj < grid_size; jj++)
        {
          double t_j = t_grid[jj];
          int status = gsl_odeiv2_driver_apply (d, &t, t_j, y);
          if (status != GSL_SUCCESS)
          {
            std::ostringstream err;
            err << "Failed to solve differential equation to compute dNeffCMB and etaBB for GeneralCosmoALP model at iteration step "<< ii <<". Invalidating point";
            invalid_point().raise(err.str());
          }
          T_evo[jj] = y[0];
        }

        gsl_odeiv2_driver_free (d);

        // update Hubble rate for next iteration
        SM.set_HT_evo(T_evo);

        eta_ratio = pow(T_evo[grid_size-1]/T_evo[0], 3) * exp(3.0*SM.get_H_int()[grid_size-1]);
        Neff_ALP = 3*pow(Tnu_grid[grid_size-1]/T_evo[grid_size-1], 4) * 3.85280407965; // (11/4)^(4/3) = 3.85280407965
        Neff_SM = 3*pow(Tnu_grid[0]/T_evo[0], 4) * 3.85280407965; // (11/4)^(4/3) = 3.85280407965
        dNeff = Neff_ALP - Neff_SM;
        ii ++;
      }

      // invalidate point if results not converged after 'max_iter' steps
      if( ((fabs(temp_eta_ratio-eta_ratio) > epsrel) || fabs(temp_dNeff - dNeff) > epsrel) && (ii >= max_iter) )
      {
        std::ostringstream err;
        err << "Computation of dNeffCMB and etaBBN for GeneralCosmoALP model did not converge after n = "<< ii <<" iterations. You can increase the maximum number of iterations with the run Option 'max_iter'. Invalidating point.";
        invalid_point().raise(err.str());
      }

      result["dNeff"] = dNeff;
      result["Neff_ratio"] = Neff_ALP/Neff_SM;
      result["eta_ratio"] = eta_ratio;
      logger() << "GeneralCosmoALP model: calculated Neff @BBN to be " << result["dNeff"] <<", and etaBB(ALP)/etaBBN(SM) = " << result["eta_ratio"] << ". Calculation converged after "<< ii <<" iterations." << EOM;

    }


  } // namespace CosmoBit

} // namespace Gambit
