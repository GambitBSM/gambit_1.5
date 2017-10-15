//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Routines to help translate between SLHA2 sfermions
///  and SLHA1 (or similar) sfermions.
///
///  Note that when family mixing occurs, there is no clear
///  definition of ~t_1, ~t_2, etc, as these are neither
///  gauge eigenstates nor mass eigenstates.  Extensive
///  matrix manipulations would be required in order to
///  define the family states as the results of diagonalising
///  a 2x2 submatrix of the full 6x6 mass matrix.  Here
///  we instead define the family states to be the two mass
///  eigenstates with the largest contributions from gauge
///  eigenstates of the appropriate family.  For example,
///  this means that ~t_1 is the lightest of the two mass
///  eigenstates that have the largest combined contribution
///  from the gauge eigenstates ~t_L and ~t_R.
///
///  *********************************************
///
///  Authors:
///
///  \author Peter Athron
///          (peter.athron@coepp.org.au)
///  \date 2015
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015
///
///  *********************************************


#ifndef __MSSM_slhahelp_hpp__
#define __MSSM_slhahelp_hpp__

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <set>

#include "gambit/Elements/subspectrum.hpp"
#include "gambit/Elements/spectrum.hpp"
#include "gambit/Utils/util_types.hpp"

namespace Gambit
{

   namespace slhahelp
   {

      /// Add a disclaimer about the absence of a MODSEL block in a generated SLHAea object
      void add_MODSEL_disclaimer(SLHAstruct& slha, const str& object);

      /// Simple helper function for adding missing SLHA1 2x2 family mixing matrices to an SLHAea object.
      void attempt_to_add_SLHA1_mixing(const str& block, SLHAstruct& slha, const str& type,
                                       const SubSpectrum& spec, double tol, str& s1, str& s2, bool pterror);

      /// ***************** Gauge <-> Mass Eigenstate Helpers ****************
      /// @{

      /// Identifies the mass eigenstate with largest gauge eigenstate content.
      ///   @{
      /// Version that tests internally agains a user-requested tolerance, either
      /// raising a GAMBIT error (if pterror_only = false) or invalidating a point.
      str mass_es_from_gauge_es(str gauge_es, const SubSpectrum& mssm, double tol,
                                str context, bool pterror_only);
      /// Version returning maximum mixing
      str mass_es_from_gauge_es(str gauge_es, double & max_mixing,
                                const SubSpectrum& mssm);
      /// Version returning gauge composition of identified mass eigenstate
      str mass_es_from_gauge_es(str gauge_es,
                                std::vector<double> & gauge_composition,
                                const SubSpectrum& mssm);
      /// Version returning maximum mixing and full gauge composition of
      /// identified mass eigenstate.
      str mass_es_from_gauge_es(str gauge_es, double & max_mixing,
                                std::vector<double> & gauge_composition,
                                const SubSpectrum& mssm);
      ///   @}

      /// Identifies the gauge eigenstate with largest mass eigenstate content.
      ///   @{
      /// Version that tests internally agains a user-requested tolerance, either
      /// raising a GAMBIT error (if pterror_only = false) or invalidating a point.
      str gauge_es_from_mass_es(str mass_es, const SubSpectrum& mssm,
                                double tol, str context, bool pterror_only);
      /// Version returning maximum mixing
      str gauge_es_from_mass_es(str mass_es, double & max_mixing,
                                const SubSpectrum& mssm);
      /// Version returning mass composition of identified gauge eigenstate
      str gauge_es_from_mass_es(str mass_es,
                                std::vector<double> & mass_composition,
                                const SubSpectrum& mssm);
      /// Version returning maximum mixing and full mass composition of
      /// identified gauge eigenstate.
      str gauge_es_from_mass_es(str mass_es, double & max_mixing,
                                std::vector<double> & mass_composition,
                                const SubSpectrum& mssm);
      ///   @}
      /// @}


      ///***************** Family <-> Mass Eigenstate Helpers ***********/
      /// @{

      /// Identifies the mass eigenstate that best matches the requested family state.
      ///   @{
      /// Version that tests internally agains a user-requested tolerance for family mixing, either
      /// raising a GAMBIT error (if pterror_only = false) or invalidating a point.
      str mass_es_closest_to_family(str familystate, const SubSpectrum& mssm,
                                    double tol, str context, bool pterror_only);
      /// Version returning mixing elements of the resulting mass eigenstate
      /// into the two gauge eigenstates of the requested family.  To test
      /// against family mixing, check that the square sum of elements of this
      /// mixing matrix row are sufficiently close to 1.  That is, compare
      /// gauge_composition(1)^2 + gauge_composition(2)^2 to 1-tolerance.
      str mass_es_closest_to_family(str familystate,
                                    std::vector<double> & gauge_composition,
                                    const SubSpectrum& mssm);
      /// Version returning the square sum of gauge mixing elements
      str mass_es_closest_to_family(str familystate,
                                    double & sum_sqr_mix,
                                    const SubSpectrum& mssm);
      /// Version returning mixing elements of the resulting mass eigenstate
      /// into the two gauge eigenstates of the requested family, and off-family mixing.
      str mass_es_closest_to_family(str familystate,
                                    std::vector<double> & gauge_composition,
                                    std::vector<double> & off_family_mixing,
                                    const SubSpectrum& mssm);
      ///   @}

      /// Identifies the family state that best matches the requested mass eigenstate.
      ///   @{
      /// Version that tests internally agains a user-requested tolerance for family mixing, either
      /// raising a GAMBIT error (if pterror_only = false) or invalidating a point.
      str family_state_closest_to_mass_es(str mass_es, const SubSpectrum& mssm,
                                          double tol, str context, bool pterror_only);
      /// Version returning the mass eigenstate composition of the gauge
      /// eigenstate that best matches the requested mass eigenstate.
      str family_state_closest_to_mass_es(str mass_es, std::vector<double> & mass_comp,
                                          const SubSpectrum& mssm);
      /// Version returning the summed squares of the contributions to the
      /// gauge eigenstate that best matches the requested mass eigenstate,
      /// of the two mass eigenstates that look most like the resulting family.
      /// (Seriously, just use the tol version.)  To test against family
      /// mixing, you can check that this square of elements is sufficiently
      /// close to 1.
      str family_state_closest_to_mass_es(str mass_es, double & sum_sqr_mix,
                                          const SubSpectrum& mssm);
      /// Version returning the mass eigenstate composition of the best-matching
      /// gauge eigenstate, and the summed squares of the contributions to this
      /// from the two mass eigenstates that look most like the resulting family.
      str family_state_closest_to_mass_es(str mass_es, double & sum_sqr_mix,
                                          std::vector<double> & mass_comp,
                                          const SubSpectrum& mssm);
      ///   @}

      /// Identifies the two mass eigenstates which best match a requested family,
      /// as well as the resulting 2x2 family mixing matrix between them.
      /// The matrix has the form (Mix_{11}, Mix_{12}, Mix_{21}, Mix_{22}).
      ///   @{
      /// Version that tests internally agains a user-requested tolerance for family mixing, either
      /// raising a GAMBIT error (if pterror_only = false) or invalidating a point.
      std::vector<double> family_state_mix_matrix(str type /*"~u", "~d" or "~e-"*/, int generation,
                                                  str & mass_es1, str & mass_es2, const SubSpectrum& mssm,
                                                  double tol, str context, bool pterror_only);
      /// Version that leaves the test up to the user.
      /// To test that there is negligible family mixing, you can check that for both rows of the
      /// family mixing matrix, the sum of squares of elements is sufficently close to 1.  That is,
      /// check Mix_{11}^2 + Mix_{12}^2 > 1-tolerance && Mix_{21}^2 + Mix_{22}^2 > 1-tolerance.
      /// where vec is the std::vector returned by this method
      std::vector<double> family_state_mix_matrix(str type /*"~u", "~d" or "~e-"*/, int generation,
                                                  str & mass_es1, str & mass_es2, const SubSpectrum& mssm);
      ///   @}

      /// @}

      /// Add an entire MSSM spectrum to an SLHAea object
      void add_MSSM_spectrum_to_SLHAea(const SubSpectrum& mssmspec, SLHAstruct& slha, int slha_version);

   }  // namespace slhahelp


   /// Structure to hold mass eigenstate pseudonyms for gauge eigenstates
   struct mass_es_pseudonyms
   {

   public:

      /// Constructor
      mass_es_pseudonyms() : filled (false) {}

      /// Known pseudonym strings
      /// @{
      ///   Particles
      ///   @{
      str isdl;
      str isul;
      str issl;
      str iscl;
      str isb1;
      str ist1;
      str isell;
      str isnel;
      str ismul;
      str isnmul;
      str istau1;
      str isntaul;
      str isdr;
      str isur;
      str issr;
      str iscr;
      str isb2;
      str ist2;
      str iselr;
      str ismur;
      str istau2;
      //str iner;   uncomment to add right-handed nus
      //str inmur;  uncomment to add right-handed nus
      //str intaur; uncomment to add right-handed nus
      ///   @}
      ///   Anti-particles
      ///   @{
      str isdlbar;
      str isulbar;
      str isslbar;
      str isclbar;
      str isb1bar;
      str ist1bar;
      str isellbar;
      str isnelbar;
      str ismulbar;
      str isnmulbar;
      str istau1bar;
      str isntaulbar;
      str isdrbar;
      str isurbar;
      str issrbar;
      str iscrbar;
      str isb2bar;
      str ist2bar;
      str iselrbar;
      str ismurbar;
      str istau2bar;
      //str inerbar;   uncomment to add right-handed nus
      //str inmurbar;  uncomment to add right-handed nus
      //str intaurbar; uncomment to add right-handed nus
      ///   @}
      /// @}

      /// Maps relating the pseudonym strings in both directions
      std::map<str,str> gauge_family_eigenstates;
      std::map<str,str> mass_eigenstates;

      /// Struct has already been filled or not.
      bool filled;

      /// Fill strings and maps in struct
      void fill(const SubSpectrum&, double, bool, bool);

      /// Refill strings and maps in struct
      void refill(const SubSpectrum&, double, bool, bool);

      /// Debug printer for pseudonyms
      void debug_print(const SubSpectrum&);

      /// Gauge state debug printer for pseudonyms
      void debug_print_gauge(const SubSpectrum&, str&, str&, double&);

      /// Family state debug printer for pseudonyms
      void debug_print_family(const SubSpectrum&, str&, str&, double&, double&);


   private:

      /// Helper functions for filling maps from MSSM gauge eigenstates to mass eigenstates.
      /// @{
      void fill_mass_es_psn_gauge(str&, str&, str, const SubSpectrum&, double, bool, bool);
      void fill_mass_es_psn_family(str&, str&, str, const SubSpectrum&, double, bool, bool);
      /// @}

    };


} // namespace gambit

#endif
