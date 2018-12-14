//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  SQLite printer class member function definitions
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

#include <iostream>
#include <sstream>
#include <sqlite3.h> // SQLite3 C interface 

// Gambit
#include "gambit/Printers/printers/sqliteprinter.hpp"

namespace Gambit
{
  namespace Printers
  {
    // Compute unique integer from two integers
    // We use this to turn MPI rank and point ID integers into an SQLite row ID
    std::size_t pairfunc(const std::size_t i, const std::size_t j)
    {
        // The Cantor pairing function should be good enough for this purpose I think
        // If we exceed the maximum size of size_t then we'll have to use a more space-efficient pairing function
        return ((i+j)*(i+j+1))/2 + j; 
    }

    // Null callback function for SQLite3 database queries where we are just inserting
    // information and not reading anything
    static int null_callback(void *NotUsed, int argc, char **argv, char **azColName) 
    {
       return 0;
    }

    // Constructor
    SQLitePrinter::SQLitePrinter(const Options& options, BasePrinter* const primary)
    : BasePrinter(primary,options.getValueOrDef<bool>(false,"auxilliary"))
    , printer_name("Primary printer")
    , database_file("uninitialised")
    , table_name("uninitialised")
    , db(NULL);
    , db_is_open(false)
    {
        // Get path of database file where results should ultimately end up
        std::ostringstream ff;
        if(options.hasKey("output_path"))
        {
          ff << options.getValue<std::string>("output_path") << "/";
        }
        else
        {
          ff << options.getValue<std::string>("default_output_path") << "/";
        }

        if(options.hasKey("output_file"))
        {
          ff << options.getValue<std::string>("output_file");
        }
        else
        {
          printer_error().raise(LOCAL_INFO, "No 'output_file' entry specified in the options section of the Printer category of the input YAML file. Please add a name there for the output sqlite database file of the scan.");
        }

        database_file = ff.str();

        // Get the name of the data table for this run
        table_name = options.getValueOrDef("table_name","results");
        
        // Create/open the database file 
        open_db(database_file);

        // Create/open the results table in the database
        open_table(table_name);
    }

    // Open database and 'attach' it to this object
    void SQLitePrinter::open_db(const std::string& path)
    { 
      // Check if we already have an open database
      if(db!=NULL)
      {
         std::stringstream err
         err << "Refused to open database file '"<<path<<"'; a database file pointer has already been attached to the SQLite printer!";
         printer_error().raise(LOCAL_INFO,err.str()); 
      }

      if(db_is_open)
      {
         std::stringstream err
         err << "Refused to open database file '"<<path<<"'; a database file is already flagged by the SQLite printer as open!";
         printer_error().raise(LOCAL_INFO,err.str());  
      }

      int rc; // return code

      rc = sqlite3_open(path, &db);
   
      if( rc ) {
         std::stringstream err
         err << "Failed to open database file '"<<path<<"':" << sqlite3_errmsg(db);
         printer_error.raise(LOCAL_INFO, err.str());
      }
      else
      {
         db_is_open = true;
      }
    }

    void SQLitePrinter::close_db()
    {
      sqlite3_close(db); // Check error code?
      db_is_open = false;
    }

    // Create results table
    void make_table(const std::string& name)
    {
      // Construct the SQLite3 statement
      std::stringstream sql;
      sql << "CREATE TABLE "<<table_name<<"("
          << "pairID   INT PRIMARY KEY NOT NULL,"
          << "MPIrank  INT             NOT NULL,"
          << "PointID  INT             NOT NULL,"
          << ");";

      /* Execute SQL statement */
      int rc;
      char *zErrMsg = 0;
      rc = sqlite3_exec(db, sql, null_callback, 0, &zErrMsg);
  
      if( rc != SQLITE_OK ){
          std::stringstream err;
          err << "SQL error: " << zErrMsg;
          printer_error.raise(LOCAL_INFO,err.str());
          sqlite3_free(zErrMsg);
      }

    }
 
  }
}
