//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  SQLite printer class declaration
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2018 Dec
///
///  *********************************************


#ifndef __sqliteprinter_hpp__
#define __sqliteprinter_hpp__

// Gambit
#include "gambit/Printers/baseprinter.hpp"

namespace Gambit
{
  namespace Printers
  {
    // Compute unique integer from two integers
    // We use this to turn MPI rank and point ID integers into an SQLite row ID
    inline std::size_t pairfunc(const std::size_t i, const std::size_t j);
 
    /// The main printer class for output to SQLite database
    class SQLitePrinter : public BasePrinter
    {
      public:
        /// Constructor (for construction via inifile options)
        HDF5Printer(const Options&, BasePrinter* const primary = NULL);

        /// Destructor
        ~SQLitePrinter() {}

        /// Virtual function overloads:
        ///@{

        // Initialisation function
        // Run by dependency resolver, which supplies the functors with a vector of VertexIDs whose requiresPrinting flags are set to true.
        void initialise(const std::vector<int>&);
        void flush();
        void reset(bool force=false);
        void finalise(bool abnormal=false);

        // Get options required to construct a reader object that can read
        // the previous output of this printer.
        Options resume_reader_options();
 
        ///@}

        ///@{ Print functions
        using BasePrinter::_print; // Tell compiler we are using some of the base class overloads of this on purpose.
        #define DECLARE_PRINT(r,data,i,elem) void _print(elem const&, const std::string&, const int, const uint, const ulong);
        BOOST_PP_SEQ_FOR_EACH_I(DECLARE_PRINT, , HDF5_TYPES)
        #ifndef SCANNER_STANDALONE
          BOOST_PP_SEQ_FOR_EACH_I(DECLARE_PRINT, , HDF5_MODULE_BACKEND_TYPES)
        #endif
        #undef DECLARE_PRINT
        ///@}

         /// Helper print functions
        // Used to reduce repetition in definitions of virtual function overloads
        // (useful since there is no automatic type conversion possible)
        template<class T>
        void template_print(T const& value, const std::string& label, const int IDcode, const unsigned int mpirank, const unsigned long pointID)
        {
 
        } 

     private:

        // Path to output SQLite database file
        std::string database_file;

        // Name of data table to store results for this run
        std::string table_name;

        // Pointer to output sqlite3 database
        sqlite3* db;

        // Bool to record if we already have a database file open
        bool db_is_open;

        // Bool to record if an output table exists yet
        bool results_table_exists;

        // Set to record whether table columns have been created 
        std::map<std::string,std::string> column_record; 

        // Check if the outbase database is open and the results table exists
        bool output_ready(); 

        // Open database and 'attach' it to this object
        // A database will be created if it doesn't exist
        void open_db(const std::string&);
    
        // Close the database file that is attached to this object
        void close_db();

        // Create results table
        void make_table(const std::string&);

        // Check that a table column exists, and create it if needed
        void ensure_column_exists(const std::string&);

        // Create table insert operation  
    }

    // Register printer so it can be constructed via inifile instructions
    // First argument is string label for inifile access, second is class from which to construct printer
    LOAD_PRINTER(sqlite, SQLitePrinter)

  }
}
#endif
