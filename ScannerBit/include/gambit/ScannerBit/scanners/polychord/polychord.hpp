//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
///  \file
///
///  ScannerBit interface to Multinest 3.10
///
///  Header file
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (ben.farmer@gmail.com)
///  \date 2015 January
///
///  *********************************************

#ifndef __polychord_hpp__
#define __polychord_hpp__

#include "gambit/ScannerBit/scanner_plugin.hpp"

// Auxilliary classes and functions needed by polychord

struct Settings
{
    int nDims;
    int nDerived;
    int nlive;
    int num_repeats;
    int nprior;
    bool do_clustering;
    int feedback;
    double precision_criterion;
    int max_ndead;
    double boost_posterior;
    bool posteriors;
    bool equals;
    bool cluster_posteriors;
    bool write_resume;
    bool write_paramnames;
    bool read_resume;
    bool write_stats;
    bool write_live;
    bool write_dead;
    bool write_prior;
    double compression_factor;
    std::string base_dir;
    std::string file_root;
    int seed;

    Settings(int _nDims=0,int _nDerived=0);
};

void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings);


namespace Gambit
{

   namespace PolyChord
   {

      /// Typedef for the ScannerBit pointer to the external loglikelihood function
      typedef Gambit::Scanner::like_ptr scanPtr;

      /// Bring printer_interface and printer names into this scope
      using Gambit::Scanner::printer_interface;
      using Gambit::Scanner::printer;

      /// Class to connect PolyChord log-likelihood function and ScannerBit likelihood function
      class LogLikeWrapper
      {
         private:
            /// Scanner pointer (points to the ScannerBit provided log-likelihood function)
            scanPtr boundLogLike;

            /// Reference to a printer_interface object
            printer_interface& boundPrinter;

            /// Number of free parameters
            int my_ndim;

            /// Variable to indicate whether the dumper function has been run at least once
            bool dumper_runonce;

         public:
            /// Constructor
            LogLikeWrapper(scanPtr, printer_interface&, int);
   
            /// Main interface function from PolyChord to ScannerBit-supplied loglikelihood function 
            double LogLike(double*, int, int);

            /// Main interface to PolyChord dumper routine   
            void dumper(int, int, int, double*, double*, double*, double, double, double);
      };

      ///@{ Plain-vanilla C-functions to pass to Multinest for the callbacks
      // Note: we are using the c interface from cwrapper.f90, so the function
      // signature is a little different than in the multinest examples.
      double callback_loglike(double*, int, int, void*);

      void callback_dumper(int, int, int, double*, double*, double*, double, double, double, void*);
      ///@}      

   } // End PolyChord namespace

} // End Gambit namespace

#endif
