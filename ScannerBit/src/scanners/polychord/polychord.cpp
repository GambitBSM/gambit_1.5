//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
///  \file
///
///  ScannerBit interface to Multinest 3.10
///
///  *********************************************
///
///  Authors (add name and date if you modify):
//
///  \author Ben Farmer
///          (ben.farmer@gmail.com)
///  \date October 2013 - Aug 2016
///
///  *********************************************

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <iomanip>  // For debugging only

#include "gambit/ScannerBit/scanner_plugin.hpp"
#include "gambit/ScannerBit/scanners/multinest/multinest.hpp"
#include "gambit/Utils/yaml_options.hpp"
#include "gambit/Utils/util_functions.hpp"


namespace Gambit
{
   namespace MultiNest
   {
      /// Global pointer to loglikelihood wrapper object, for use in the MultiNest callback functions
      LogLikeWrapper *global_loglike_object;
   }
}

/// Typedef for the ScannerBit pointer to the external loglikelihood function
typedef Gambit::Scanner::like_ptr scanPtr;


/// =================================================
/// Interface to ScannerBit
/// =================================================

scanner_plugin(multinest, version(3, 10))
{
   // An error is thrown if any of the following entries are not present in the inifile (none absolutely required for MultiNest).
   reqd_inifile_entries();

   // Tell cmake system to search known paths for these libraries; any not found must be specified in config/scanner_locations.yaml.
   reqd_libraries("nest3");

   // Pointer to the (log)likelihood function
   scanPtr LogLike;

   /// The constructor to run when the MultiNest plugin is loaded.
   plugin_constructor
   {
      // Retrieve the external likelihood calculator
      LogLike = get_purpose(get_inifile_value<std::string>("like"));
      if (LogLike->getRank() == 0) std::cout << "Loading MultiNest nested sampling plugin for ScannerBit." << std::endl;
   }

   /// The main routine to run for the MultiNest scanner.
   int plugin_main (void)
   {
      /// ************
      /// TODO: Replace with some wrapper? Maybe not, this is already pretty straightforward,
      /// though perhaps a little counterintuitive that the printer is the place to get this
      /// information.
      bool resume_mode = get_printer().resume_mode();
      /// ************

      // Retrieve the dimensionality of the scan.
      int ma = get_dimension();

      // Retrieve the global option specifying the minimum interesting likelihood
      double gl0 = get_inifile_value<double>("likelihood: model_invalid_for_lnlike_below");
      // Retrieve the global option specifying the likelihood offset to use
      double offset = get_inifile_value<double>("likelihood: lnlike_offset", 0.);
      // Make sure the likleihood functor knows to apply the offset internally in ScannerBit
      LogLike->setPurposeOffset(offset);
      // Offset the minimum interesting likelihood by the offset
      gl0 = gl0 + offset;


      // PolyChord algorithm options.
      Settings settings;
      settings.nDims = ma;
      settings.nDerived = 0;
      settings.nlive = get_inifile_value<int>("nlive", nDims*25);                  // number of live points
      settings.num_repeats = get_inifile_value<int>("num_repeats", nDims*5);       // length of slice sampling chain
      settings.nprior = get_inifile_value<int>("nprior", nlive*10);                // number of prior samples to begin algorithm with
      settings.do_clustering = get_inifile_value<bool>("do_clustering", true);     // Whether or not to perform clustering
      settings.feedback = get_inifile_value<int>("feedback", 1);                   // Feedback level
      settings.precision_criterion = get_inifile_value<double>("tol", 0.5);      // Stopping criterion (consistent with multinest)
      settings.max_ndead = get_inifile_value<double>("maxiter", 0);                  // Max no. of iterations, a non-positive value means infinity (consistent with multinest).
      settings.boost_posterior = get_inifile_value<double>("boost_posterior",0.); // Increase the number of posterior samples produced
      bool outfile (get_inifile_value<bool>("outfile", true));                // write output files?
      settings.posteriors = outfile;
      settings.equals = outfile;
      settings.cluster_posteriors = outfile;
      settings.write_paramnames = outfile;
      settings.write_stats = outfile;
      settings.write_live = outfile; 
      settings.write_dead = outfile;
      settings.write_prior = outfile;
      settings.write_resume = resume_mode;
      settings.read_resume = resume_mode;
      settings.compression_factor = get_inifile_value<double>("compression_factor",0.36787944117144233);
      settings.base_dir = Gambit::Utils::ensure_path_exists(get_inifile_value<std::string>("default_output_path")+"PolyChord");
      settings.file_root = get_inifile_value<std::string>("root", "native");

      settings.seed = get_inifile_value<int>("seed",-1);


      if(resume==1 and outfile==0)
      {
        // It is stupid to be in resume mode while not writing output files.
        // Means subsequent resumes will be impossible. Throw an error.
        scan_error().raise(LOCAL_INFO,"Error from PolyChord ScannerBit plugin! Resume mode is activated, however "
                                      "PolyChord native output files are set to not be written. These are needed "
                                      "for resuming; please change this setting in your yaml file (set option \"outfile: 1\")");
      }

      // Setup auxilliary streams. These are only needed by the master process,
      // so let's create them only for that process
      int myrank = get_printer().get_stream()->getRank(); // MPI rank of this process
      if(myrank==0)
      {
         // Get inifile options for each print stream
         Gambit::Options txt_options   = get_inifile_node("aux_printer_txt_options");
         //Gambit::Options stats_options = get_inifile_node("aux_printer_stats_options"); //FIXME
         Gambit::Options live_options  = get_inifile_node("aux_printer_live_options");

         // Options to desynchronise print streams from the main Gambit iterations. This allows for random access writing, or writing of global scan data.
         //stats_options.setValue("synchronised",false); //FIXME
         txt_options.setValue("synchronised",false);
         live_options.setValue("synchronised",false);

         // Initialise auxiliary print streams
         get_printer().new_stream("txt",txt_options);
         //get_printer().new_stream("stats",stats_options); //FIXME
         get_printer().new_stream("live",live_options);
      }

      // Ensure that MPI processes have the same IDs for auxiliary print streams;
      Gambit::Scanner::assign_aux_numbers("Posterior","LastLive");

      // Create the object that interfaces to the MultiNest LogLike callback function
      Gambit::MultiNest::LogLikeWrapper loglwrapper(LogLike, get_printer(), ndims);
      Gambit::MultiNest::global_loglike_object = &loglwrapper;

      //Run MultiNest, passing callback functions for the loglike and dumper.
      if(myrank == 0) std::cout << "Starting MultiNest run..." << std::endl;
      polychord_c_interface( 
              Gambit::PolyChord::callback_loglike, 
              Gambit::PolyChord::callback_prior, 
              nlive, 
              num_repeats,
              nprior,
              do_clustering,
              feedback,
              precision_criterion,
              max_ndead,
              boost_posterior,
              posteriors,
              equals,
              cluster_posteriors,
              write_resume,
              write_paramnames,
              read_resume,
              write_stats,
              write_live,
              write_dead,
              write_prior,
              compression_factor, 
              nDims,
              nDerived,
              base_dir,
              file_root,
              nGrade,
              grade_frac,
              grade_dims,
              n_nlives,
              loglikes,
              nlives,
              seed
                  );



      polychord_c_interface(
              Gambit::PolyChord::callback_loglike, Gambit::PolyChord::callback_prior, IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol,
              root, seed, pWrap, fb, resume, outfile, initMPI, ln0, maxiter,
              Gambit::MultiNest::callback_loglike, Gambit::MultiNest::callback_dumper, context);
      if(myrank == 0) std::cout << "Multinest run finished!" << std::endl;
      return 0;

   }

}


/// =================================================
/// Function definitions
/// =================================================

namespace Gambit {

   namespace MultiNest {

      ///@{ Plain-vanilla functions to pass to PolyChord for the callback
      // Note: we are using the c interface from cwrapper.f90, so the function
      // signature is a little different than in the multinest examples.
      double callback_loglike(double *Cube, int ndim, double* phi, int nderived)
      {
         // Call global interface to ScannerBit loglikelihood function
         // Could also pass this object in via context pointer, but that
         // involves some casting and could risk a segfault.
         return global_loglike_object->LogLike(Cube, ndim, phi nderived);
      }

      void callback_dumper(int nSamples, int nlive, int nPar, double *physLive,
                           double *posterior, double *paramConstr,
                           double maxLogLike, double logZ, double logZerr,
                           void*)
      {
         global_loglike_object->
            dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr,
                   maxLogLike, logZ, logZerr);
      }
      ///@}


      /// LogLikeWrapper Constructor
      LogLikeWrapper::LogLikeWrapper(scanPtr loglike, printer_interface& printer, int ndim)
        : boundLogLike(loglike), boundPrinter(printer), my_ndim(ndim), dumper_runonce(false)
      { }

      /// Main interface function from PolyChord to ScannerBit-supplied loglikelihood function
      /// This is the function that will be passed to PolyChord as the
      /// loglike callback routine
      ///
      /// Input arguments
      /// ndim    = dimensionality (total number of free parameters) of the problem
      /// npars   = total number of free plus derived parameters
      ///
      /// Input/Output arguments
      /// Cube[npars]  = on entry has the ndim parameters in unit-hypercube
      ///                on exit, the physical parameters plus copy any derived parameters
      ///                you want to store with the free parameters
      ///
      /// Output arguments
      /// lnew = loglikelihood
      ///
      double LogLikeWrapper::LogLike(double *Cube, int ndim, int)
      {
         //convert C style array to C++ vector class
         std::vector<double> unitpars(Cube, Cube + ndim);

         double lnew = boundLogLike(unitpars);

         // Done! (lnew will be used by PolyChord to guide the search)
         return lnew;
      }
   }

}

