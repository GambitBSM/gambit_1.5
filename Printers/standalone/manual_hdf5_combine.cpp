//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  This is code for a stand-alone binary which
///  performs the end-of-run combination of
///  temporary HDF5 files from a GAMBIT scan.
///
///  It is essentially the same code that GAMBIT
///  runs automatically in the HDF5 printer, just
///  moved into a stand-alone tool so that the
///  combination can be done "manually" if needed.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2018 Mar
///
///  *********************************************

#include <stdlib.h>
#include <getopt.h>
#include <vector>
#include <sstream>
#include <utility>
#include <string>
#include <sys/stat.h>

// GAMBIT headers
#include "gambit/Printers/printers/hdf5printer/hdf5_combine_tools.hpp"
#include "gambit/Utils/stream_overloads.hpp"
#include "gambit/Utils/util_functions.hpp"

// Annoying other things we need due to mostly unwanted dependencies
#include "gambit/Utils/static_members.hpp"

using namespace Gambit;
using namespace Printers;

void usage()
{
    std::cout << "\nusage: hdf5combine <filename> <group> [options] "
          "\n "
          "\n  filename - Base name of the temporary files to be combined" 
          "\n               e.g. /path/to/my/hdf5file.hdf5 "
          "\n  group    - location inside the HDF5 files which contains the datasets to be combined"
          "\n"
          "\nOptions:"
          "\n   -h/--help             Display this usage information"
          "\n   -f/--force            Attempt combination while ignoring missing temporary files"
          "\n   -c/--cleanup          Delete temporary files after successful combination"
          "\n   -o/--out <path>       Set output folder (default is same folder as temp files)"
          "\n                         (note, this first creates the output in the default place,"
          "\n                         and then just moves it)"
          "\n\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
   // Process command line options
   // Must at least have two arguments.
   if (argc < 3)
   {
      usage();
   }

   // First, the required options (folder and temporary file base name)
   std::string finalfile  = argv[1];
   std::string group      = argv[2];
   std::string outpath;

   const struct option primary_options[] =
   {
     {"force", no_argument, 0, 'f'},
     {"help",  no_argument, 0, 'h'},
     {"cleanup", no_argument, 0 , 'c'},
     {"out", required_argument, 0, 'o'},
     {0,0,0,0},
   };

   int index;
   int iarg=0;
   std::string filename;

   bool error_if_inconsistent = true;
   bool do_cleanup = false;
   bool move_output = false;

   while(iarg != -1)
   {
     iarg = getopt_long(argc, argv, "fhco:", primary_options, &index);
     switch (iarg)
     {
       case 'h':
         // Display usage information and exit
         usage();
         break;
       case 'f':
         // Ignore missing temporary files
         error_if_inconsistent = false;
         break;
       case 'c':
         do_cleanup = true;
         break;
       case 'o':
         outpath = optarg; 
         move_output = true;
         break;
       case '?':
         // display usage message and quit (also happens on unrecognised options)
         usage();
         break;
     }
   }

   // Report what we are going to do
   std::cout <<"  Searching for temporary files to combine" <<std::endl;
   std::cout <<"     with base filename: "<<finalfile<<std::endl;



   // Search for the temporary files to be combined
   std::pair<std::vector<std::string>,std::vector<int>> out = HDF5::find_temporary_files(finalfile);
   std::vector<std::string> tmp_files = out.first;
   std::vector<int> missing = out.second;

   std::cout <<"  Found "<<tmp_files.size()<<" temporary files to be combined"<<std::endl;
   // Check if all temporary files found (i.e. check if output from some rank is missing)
   if(missing.size()>0)
   {
     std::cout <<"  WARNING! Temporary files missing from the following ranks: "<<missing<<std::endl; 
     if(error_if_inconsistent)
     {
       std::ostringstream errmsg;
       errmsg << "  Could not locate all the expected temporary output files (found "<<tmp_files.size()<<" temporary files, but are missing the files from the following ranks: "<<missing<<")! If you want to attempt the combination anyway, please add the -f flag to the command line invocation.";
       std::cerr << errmsg.str() << std::endl;
       exit(EXIT_FAILURE); 
     }
   }

   // Name of temporary combined file, if one exists
   std::ostringstream name;
   name << finalfile << "_temp_combined";
   std::string tmp_comb_file = name.str();

   bool combined_file_exists = Utils::file_exists(tmp_comb_file);
   if(combined_file_exists)
   {
      std::cout<<"  Found pre-existing (temporary) combined file"<<std::endl;
   } 
   else 
   {
      std::cout<<"  No pre-existing (temporary) combined file found"<<std::endl;  
   }

   int num = tmp_files.size(); // We don't actually use their names here, Greg's code assumes that they
                               // follow a fixed format and they all exist. We check for this before
                               // running this function, so this should be fine.

   HDF5::combine_hdf5_files(tmp_comb_file, finalfile, group, num, combined_file_exists, do_cleanup);

   std::cout<<"  Combination finished successfully!"<<std::endl;
   if(move_output)
   {
      std::string fname(Utils::base_name(tmp_comb_file));
      std::string outname(outpath+"/"+fname);
      std::cout<<"  Moving output from "<<std::endl
               <<"  "<<tmp_comb_file<<std::endl
               <<"  to"<<std::endl
               <<"  "<<outname<<std::endl;
      int status = std::system(("mv " + tmp_comb_file+ " " + outname).c_str());
      if(WEXITSTATUS(status)!=0)
      {
         std::cerr << "  FAILED to move output to requested directory!! It should still be located at"<<std::endl
                   << "    "<<tmp_comb_file<<std::endl;
      }
   }

   // Finished!
   exit(EXIT_SUCCESS);   
}

