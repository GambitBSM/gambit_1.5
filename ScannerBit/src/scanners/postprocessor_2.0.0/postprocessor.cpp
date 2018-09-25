//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
///  \file
///
///  "Postprocessing" scanner plugin. Reads points
///  from old scan output and re-runs a likelihood
///  containing plugin for those same point.
///  Can perform some simple addition/subtraction
///  operations of likelihood components from
///  the new plugin output.
///
///  This is version 2 of the postproccessor; it
///  distributes the postprocessing workload by
///  a completely different algorithm to version 1.
///  This version employs a master/worker model,
///  with the master processes distributing points
///  in batches to the worker processes on request.
///  Batch size can be set via options (use batch
///  size of 1 for very slow likelihoods, use
///  large batch size for very fast likelihoods).
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2018, Sep
///
///  *********************************************

// STL
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

// GAMBIT
#include "gambit/Utils/mpiwrapper.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/new_mpi_datatypes.hpp"
#include "gambit/ScannerBit/scanners/postprocessor_2.0.0/postprocessor.hpp"
#include "gambit/ScannerBit/objective_plugin.hpp"
#include "gambit/ScannerBit/scanner_plugin.hpp"

using namespace Gambit;
using namespace Gambit::PostProcessor;

// Forward declare this template specialisation as extern so that we use the definition compiled into baseprinter.cpp
extern template std::size_t Gambit::Printers::getTypeID<double>();

// The reweigher Scanner plugin
scanner_plugin(postprocessor, version(2, 0, 0))
{
  reqd_inifile_entries("like","reader");

  /// The likelihood container plugin
  like_ptr LogLike;

  /// MPI data
  int numtasks;
  int rank;

  // Tags for messages
  const int request_work_tag=10;

  /// The reader object in use for the scan
  Gambit::Printers::BaseBaseReader* reader;

  /// The main postprocessing driver object
  PPDriver driver;

  /// Options for PPDriver;
  PPOptions settings;

  // Retrieve an integer from an environment variable
  int getintenv(const std::string& name)
  {
     int x;
     if(const char* env_p = std::getenv(name.c_str()))
     {
       std::stringstream env_s(env_p);
       env_s >> x;
       if (!env_s)
       {
          std::ostringstream err;
          err << "Tried to retrieve value of environment variable "<<name<<" as an integer, but conversion failed! String retrieved was '"<<env_s.str()<<"'";
          scan_error().raise(LOCAL_INFO,err.str());
       }
     }
     else
     {
       std::ostringstream err;
       err << "Tried to retrieve value of environment variable "<<name<<" as an integer, but it does not seem to be defined!";
       scan_error().raise(LOCAL_INFO,err.str());
     }
     return x;
  }

  /// The constructor to run when the plugin is loaded.
  plugin_constructor
  {
    int s_numtasks;
    int s_rank;

    // Get MPI data. No communication is needed, we just need to know how to
    // split up the workload. Just a straight division among all processes is
    // used, nothing fancy.
#ifdef WITH_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &s_numtasks); // MPI requires unsigned ints here, so we'll just convert afterwards
    MPI_Comm_rank(MPI_COMM_WORLD, &s_rank);

#else
    s_numtasks = 1;
    s_rank = 0;
#endif
    numtasks = s_numtasks;
    rank = s_rank;

    if(rank==0) std::cout << "Initialising 'postprocessor' plugin for ScannerBit..." << std::endl;

    // Get options for setting up the reader (these live in the inifile under:
    // Scanners:
    //  scanners:
    //    scannername:
    //      reader
    Gambit::Options reader_options = get_inifile_node("reader");

    // Initialise reader object
    get_printer().new_reader("old_points",reader_options);

    // Retrieve the reader object
    reader = get_printer().get_reader("old_points");

    // Get names of all the output data labels
    settings.data_labels = reader->get_all_labels();

    // Set up other options for the plugin
    settings.update_interval = get_inifile_value<std::size_t>("update_interval", 1000);
    settings.add_to_logl = get_inifile_value<std::vector<std::string>>("add_to_like", std::vector<std::string>());
    settings.subtract_from_logl = get_inifile_value<std::vector<std::string>>("subtract_from_like", std::vector<std::string>());
    settings.reweighted_loglike_name = get_inifile_value<std::string>("reweighted_like");

    settings.renaming_scheme = get_inifile_value<std::map<std::string,std::string>>("rename",
                          std::map<std::string,std::string>());

    settings.cut_less_than = get_inifile_value<std::map<std::string,double>>("cut_less_than",
                          std::map<std::string,double>());

    settings.cut_greater_than = get_inifile_value<std::map<std::string,double>>("cut_greater_than",
                          std::map<std::string,double>());

    settings.discard_points_outside_cuts = get_inifile_value<bool>("discard_points_outside_cuts", false);

    // Use virtual rank system?
    if(get_inifile_value<bool>("use_virtual_rank",false))
    {
        #ifdef WITH_MPI
        if(numtasks>1)
        {
          std::ostringstream err;
          err << "You have set the 'use_virtual_rank' option for the postprocessor scanner plugin to 'true', which will allow the plugin to act as if it is part of an MPI ensemble when it really isn't, however you are also running this task in an MPI batch with size > 1! You cannot use the virtual rank system at the same time as running a real MPI job! Please choose one configuration or the other and rerun the job.";
          scan_error().raise(LOCAL_INFO,err.str());
        }
        #endif
        rank     = getintenv("RANK");
        numtasks = getintenv("SIZE");
        if(rank>=numtasks)
        {
          std::ostringstream err;
          err << "Environment variable RANK was larger than permitted by SIZE ("<<numtasks<<">="<<rank<<") while running postprocessor scanner plugin with 'use_virtual_rank=true' option. This is not a valid MPI configuration, so it is an illegal choice of virtual configuration.";
          scan_error().raise(LOCAL_INFO,err.str());
        }
    }
    // Transfer MPI variables to PPOptions
    settings.rank = rank;
    settings.numtasks = numtasks;

    // Size of chunks to be distributed to worker processes
    settings.chunksize = get_inifile_value<std::size_t>("batch_size",1);

    // Finally, there is the 'Purpose' value of the likelihood container. This may well clash
    // with the old name used in the input file, so better check for this and make the user
    // change their choice if so.
    settings.logl_purpose_name = get_inifile_value<std::string>("like");
    settings.discard_old_logl = get_inifile_value<bool>("permit_discard_old_likes",false);

    // Retrieve the external likelihood calculator
    LogLike = get_purpose(settings.logl_purpose_name);

    // Do not allow GAMBIT's own likelihood calculator to directly shut down the scan.
    // This scanner plugin will assume responsibility for this process, triggered externally by
    // the 'plugin_info.early_shutdown_in_progress()' function.
    LogLike->disable_external_shutdown();

    // Path to save resume files
    std::string defpath = get_inifile_value<std::string>("default_output_path");
    settings.root = Utils::ensure_path_exists(defpath+"/postprocessor/resume");
    if(rank==0) std::cout << "root: " << settings.root << std::endl;
  }

  /// Main run function
  int plugin_main()
  {
    if(rank==0) std::cout << "Running 'postprocessor' plugin for ScannerBit." << std::endl;

    // Set up our MPI communicator
    #ifdef WITH_MPI
    GMPI::Comm ppComm;
    ppComm.dup(MPI_COMM_WORLD,"PostprocessorComm"); // duplicates MPI_COMM_WORLD
    // Message tag definitons in PPDriver class:
    #endif

    /// Determine what data needs to be copied from the input file to the new output dataset
    // Get labels of functors listed for printing from the primary printer.
    settings.all_params = get_printer().get_stream()->getPrintList();
    // There are some extra items that will also be automatically printed in all scans,
    // so we need to avoid copying those:
    settings.all_params.insert("unitCubeParameters"); // It would be better to keep the originals here, but currently cannot turn off the printing from within like_ptr.
    settings.all_params.insert("MPIrank"); // These should be re-printed the same as they were anyway
    settings.all_params.insert("pointID");
    settings.all_params.insert(settings.logl_purpose_name); // If there is a name clash and the run was not aborted, we are to discard the old data under this name.
    settings.all_params.insert(settings.reweighted_loglike_name); //   "  "
    #ifdef WITH_MPI
    settings.comm = &ppComm;
    #endif

    // Construct the main driver object
    driver = PPDriver(reader,get_printer().get_stream(),LogLike,settings);

    // Check that the supplied settings make sense
    driver.check_settings();

    // Points which have already been processed in a previous (aborted) run
    ChunkSet done_chunks; // Empty by default

    // Ask the printer if this is a resumed run or not, and check that the necessary files exist if so.
    bool resume = get_printer().resume_mode();
 
    // Vector to record which processes have been told by the master to stop. 
    // Master cannot stop until all other processes have stopped.
    std::vector<bool> process_has_stopped(numtasks);

    // Rank 0 needs to figure out which points are already processesed (if resuming)
    if(resume)
    { 
        if(rank==0)
        {
            std::cout << "Analysing previous output to determine remaining postprocessing work (may take a little time for large datasets)..." << std::endl;
            done_chunks = get_done_points(settings.root);
            std::cout << "Distributing information about remaining work to all processes..." << std::endl;
        }

        // Need to distribute these to all processes
        // It is a bit hard to distribute them in one message, so we will do it
        // one chunk at a time. Hopefully this isn't a big deal in terms of the
        // delivery times. TODO: review this if startup is too slow.
        for(ChunkSet::iterator chunk=done_chunks.begin();
               chunk!=done_chunks.end(); ++chunk)
        {
            std::size_t chunkdata[3]; // Raw form of chunk information
            chunkdata[0] = chunk->start;
            chunkdata[1] = chunk->end;
            chunkdata[2] = chunk->eff_length;
   
            ppComm.Bcast(chunkdata, 0); // Broadcast to all workers from master

            if(rank!=0)
            {
               Chunk newchunk;
               newchunk.start      = chunkdata[0];
               newchunk.end        = chunkdata[1];
               newchunk.eff_length = chunkdata[2];
               done_chunks.insert(newchunk);
            }
        }

        if(rank==0) 
        {
            std::cout << "Postprocessing resume analysis completed." << std::endl;
        }
    }

    // Tell the driver routine what points it can automatically skip
    driver.set_done_chunks(done_chunks);

    //MAIN LOOP HERE
    bool continue_processing = true;
    #ifdef WITH_MPI
      bool quit_flag_seen = false;
    #endif
    while(continue_processing)
    {
       Chunk mychunk; // Work to be performed this loop

       #ifdef WITH_MPI
         bool I_am_finished = false;
         if(rank==0 and numtasks==1)
         {
            // Compute new work for this one process.
            mychunk = driver.get_new_chunk(); //TODO: Rewrite this function. Make sure to update done_chunks as we go.
         }
         else if(rank==0)
         { 
            // Master checks for work requests from other processes
            for(int worker=1; worker<numtasks; worker++)
            {
               bool needs_work = ppComm.Iprobe(worker, request_work_tag);
               if(needs_work)
               {
                  // Receive the work request message (no information, just cleaning up)
                  int nullint = 0; // Buffer for null message
                  ppComm.Recv(&nullint,1,worker,request_work_tag);

                  // Compute new work assignment
                  Chunk newchunk = driver.get_new_chunk(); //TODO: Rewrite this function.
                  // TODO: Make sure to return Chunk(0,0) to trigger end of processing.

                  // Send work assignment
                  std::size_t chunkdata[3]; // Raw form of chunk information
                  chunkdata[0] = newchunk.start;
                  chunkdata[1] = newchunk.end;
                  chunkdata[2] = newchunk.eff_length;
                  ppComm.Send(&chunkdata, 3, worker, request_work_tag);

                  // Check if we just sent the 'stop' signal.
                  if(newchunk==stopchunk)
                  {
                     process_has_stopped[worker] = true;
                  }
               }
            }

            // Set zero-length chunk for master
            bool any_still_running=false;
            for(std::size_t i=1; i<numtasks; i++)
            {
               if(process_has_stopped[i]==false) any_still_running=true;
            }

            if(any_still_running)
            {
               mychunk = Chunk(1,1); // Zero-length chunk; master doesn't process anything, but need to continue looping
            }
            else
            {
               // Everyone has been told to stop! So now master should stop too.
               mychunk = stopchunk;
            }
         }
         else
         {
            // Worker processes request more work from master
            int nullint = 0; // Buffer for null message
            ppComm.Send(&nullint,1,0,request_work_tag);

            // Receive the work assignment
            std::size_t chunkdata[3]; // Raw form of chunk information
            ppComm.Recv(&chunkdata,3,0,request_work_tag);

            // Check if any work in the work assignment
            // If start and end are both zero then take this as the signal that
            // we are finished
            mychunk.start      = chunkdata[0];
            mychunk.end        = chunkdata[1];
            mychunk.eff_length = chunkdata[2];
         }
       #else
         // Compute new work for this one process.
         mychunk = driver.get_new_chunk(); //TODO: Rewrite this function. Make sure to update done_chunks as we go.
       #endif
       std::cout << "Rank "<<rank<<": Chunk to process is ["<<mychunk.start<<", "<<mychunk.end<<"; eff_len="<<mychunk.eff_length<<"]"<<std::endl;

       if(mychunk==stopchunk)
       {
          // Finished!
          continue_processing = false;
       }

       if(mychunk.start > mychunk.end)
       {
          // End after start, error!
          std::ostringstream err;
          err << "Work assignment for rank "<<rank
              <<" process is invalid! Chunk end ("<<mychunk.end
              <<") is before the chunk start ("<<mychunk.start
              <<")! ended due to encountering the end of the input file."
              <<" This indicates a bug in the postprocessor (or some "
              <<"bizarre corruption of the MPI message). Please report"
              <<"this.";
          std::cerr << err.str() << std::endl;
          scan_error().raise(LOCAL_INFO,err.str());
       }

       int exit_code;
       if(continue_processing)
       {
          // 0 - Finished processing all the points we were assigned
          // 1 - Saw quit flag and so stopped prematurely
          // 2 - Encountered end of input file unexpectedly
          exit_code = driver.run_main_loop(mychunk);
          std::cout << "Rank "<<rank<<": exited loop with code "<<exit_code<<std::endl;
       }
       else
       {
          // No points assigned, in shutdown mode
          exit_code=0;
       }

       if(exit_code==0)
       {
          #ifdef WITH_MPI
            I_am_finished = true;
          #endif
          std::cout << "Rank "<<rank<<" has finished processing its batch." << std::endl;
       }
       else if(exit_code==1)
       {
          // Saw quit flag, time to stop
          #ifdef WITH_MPI
            quit_flag_seen = true;
          #endif
          continue_processing = false;
       }
       else if(exit_code==2)
       {
          // That shouldn't happen; warning
          std::ostringstream err;
          err << "Postprocessing on "<<rank<<" ended due to encountering the end of the input file. This indicates that it was told to process more points than existed in the input file, which indicates a bug in the postprocessor. Your output may still be fine, but please report this bug.";
          std::cerr << err.str() << std::endl;
          scan_error().raise(LOCAL_INFO,err.str());
       }
       else
       {
          std::ostringstream err;
          err << "Postprocessing on "<<rank<<" terminated with an unrecognised return code ("<<exit_code<<"). This indicates a bug in the postprocessor, please report it.";
          scan_error().raise(LOCAL_INFO,err.str());
       }
       
    }
    //if(rank==0) std::cout << "Done!" << std::endl;
    std::cout << "Rank "<< rank<< ": Done!" << std::endl;

    // Test barrier to see if everyone made it
    #ifdef WITH_MPI
    ppComm.Barrier();
    if(rank==0) std::cout << "Passed final PP barrier" << std::endl;
    #endif
    return 0;
  }
}
