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

    // Return a comma if iterator is the last in the supplied iterable
    template <typename Iter, typename Cont>
    std::string comma_unless_last(Iter it, Cont c)
    { 
       std::string out("");
       if(it + 1 != c.end()) out = ",";
       return out;
    }

    /* Callback function for retrieving the column names and types of a table
     * Called for each row of the results table
     *
     * Arguments:
     *
     *     list - Pointer to a map from column names to types
     *    count - The number of columns in the result set
     *     data - The row's data
     *  columns - The column names
     */
    static int col_name_callback(void *colmap_in, int count, char **data, char **columns)
    {
        typedef std::map<std::string, std::string> mymaptype;
        mymaptype *colmap = static_cast<mymaptype*>(colmap_in);
   
        // We know that the column name is the second column of the results set, and the
        // type is the third column
        std::string column_name(data[1]);
        std::string column_type(data[2]);

        // Add to map
        (*colmap)[column_name] = column_type;
    
        return 0;
    }

    // Constructor
    SQLitePrinter::SQLitePrinter(const Options& options, BasePrinter* const primary)
    : BasePrinter(primary,options.getValueOrDef<bool>(false,"auxilliary"))
    , database_file("uninitialised")
    , table_name("uninitialised")
    , db(NULL)
    , db_is_open(false)
    , results_table_exists(false)
    , column_record()
    , max_buffer_length(options.getValueOrDef<std::size_t>(1,"buffer_length"))
    , buffer_info()
    , buffer_header() 
    , transaction_data_buffer()
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

    void SQLitePrinter::initialise(const std::vector<int>&)
    {
        // Don't need to initialise anything for this printer
    }

    void SQLitePrinter::reset(bool force=false)
    {
        // Not sure if this is really allowed anymore in GAMBIT.
        // Let's throw an error if it gets called.
        std::stringstream err
        err<<"reset() function in SQLitePrinter was called! This printer is not currently capable of resetting the output!";
        printer_error().raise(LOCAL_INFO,err.str()); 
    }

    void SQLitePrinter::finalise(bool /*abnormal*/);
    {
        // Dump buffer to disk. Nothing special needed for early shutdown.
        SQLitePrinter::dump_buffer();
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
   
        if( rc ) 
        {
            std::stringstream err
            err << "Failed to open database file '"<<path<<"':" << sqlite3_errmsg(db);
            printer_error().raise(LOCAL_INFO, err.str());
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

    // Make sure the outbase database is open and the results table exists
    // Throws an error if it isn't
    void SQLitePrinter::require_output_ready()
    {
        if(db==NULL || !db_is_open || !results_table_exists)
        {
            std::stringstream err
            // Something was not ready, check further and throw an error
            if(db==NULL)
            {
                err << "Output readiness check failed! Database pointer was NULL!";
                printer_error().raise(LOCAL_INFO,err.str());  
            }   

            if(!db_is_open)
            {
                err << "Output readiness check failed! Database is not flagged as open!";
                printer_error().raise(LOCAL_INFO,err.str());  
            }

            if(!results_table_exists)
            {
                err << "Output readiness check failed! Results table is not flagged as existing!";
                printer_error().raise(LOCAL_INFO,err.str());  
            }
        }
        // Else we are good to go!
    } 

    // Create results table
    void SQLitePrinter::make_table(const std::string& name)
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
        rc = sqlite3_exec(db, sql.str().c_str(), NULL, NULL, &zErrMsg);
  
        if( rc != SQLITE_OK ){
            std::stringstream err;
            err << "SQL error: " << zErrMsg;
            sqlite3_free(zErrMsg);
            printer_error().raise(LOCAL_INFO,err.str());
       }
    }

    // Check that a table column exists with the correct type, and create it if needed
    void SQLitePrinter::ensure_column_exists(const std::string& sql_col_name, const std::string& sql_col_type)
    {
        require_output_ready();
        auto it = column_record.find(sql_col_name);
        if(it == column_record.end())
        {
            // Column not marked as existing. But it might have been
            // created by another process, so we need to check the
            // database directly. It seems like the best way to do
            // this is to just attempt to add the column. If it fails
            // we can then explicitly check the column names to
            // make sure that the reason for failure was because
            // the column already existed, and not some other reason.

            std::stringstream sql;
            sql<<"ALTER TABLE "<<table_name<<" ADD COLUMN "<<sql_col_name<<" "<<sql_col_type<<";";

            /* Execute SQL statement */
            int rc;
            char *zErrMsg = 0;
            rc = sqlite3_exec(db, sql.str().c_str(), NULL, NULL, &zErrMsg);
  
            if( rc != SQLITE_OK ){
                // Operation failed for some reason. Probably because the column already
                // exists, but we better make sure.

                std::stringstream sql2;
                sql2<<"PRAGMA table_info("<<table_name<<");";

                /* Execute SQL statement */
                int rc2;
                char *zErrMsg2 = 0;
                std::map<std::string, std::string> colnames; // Will be passed to and filled by the callback function
                rc2 = sqlite3_exec(db, sql2.str().c_str(), col_name_callback, &colnames, &zErrMsg2);
 
                if( rc2 != SQLITE_OK ){
                    std::stringstream err;
                    err << "Failed to check SQL column names in output table, after failing to add a new column '"<<sql_col_name<<"' to that table."<<std::endl; 
                    err << "  First  SQL error was: " << zErrMsg << std::endl;
                    err << "  Second SQL error was: " << zErrMsg2 << std::endl;
                    sqlite3_free(zErrMsg);
                    sqlite3_free(zErrMsg2);
                    printer_error().raise(LOCAL_INFO,err.str());
                }

                // Operation successful, check if our column name exists and has the correct type
                auto jt = colnames.find(sql_col_name);
                if(jt==colnames.end())
                {
                    // Column not found
                    std::stringstream err;
                    err << "Failed to add new column '"<<sql_col_name<<"' to output SQL table! The ALTER TABLE operation failed, however it was not because the column already existed (we successfully checked and the column was not found). The SQL error was: " << zErrMsg << std::endl; 
                    sqlite3_free(zErrMsg);
                    printer_error().raise(LOCAL_INFO,err.str());
                }
                else if(jt->second != sql_col_type);
                {
                    // Column found, but has the wrong type
                    std::stringstream err;
                    err << "Failed to add new column '"<<sql_col_name<<"' to output SQL table! The column already exists, but it has the wrong type (existing column has type '"<<jt->second<<"', but we expected it to have type '"<<sql_col_type<<"'!";
                    sqlite3_free(zErrMsg);
                    printer_error().raise(LOCAL_INFO.err.str());
                }

                // Column exists and has the right type! So everything is ok after all.
            }

            // Column should exist now. Need to add the fact of this columns existence to our internal record.
            column_record[sql_col_name] = sql_col_type;     
        } 
        else if(it->second != sql_col_type)
        {
            // Records say column exists, but not with the type requested!
            std::stringstream err;
            err << "SQLitePrinter records indicated that the column '"<<sql_col_name<<"' already exists in the output table, but with a different type than has been requested (existing type is '"<<it->second<<"', requested type was '"<<sql_col_type<<"'). This indicates either duplicate names in the printer output, or an inconsistency in how the print commands have been issued.";
            printer_error().raise(LOCAL_INFO.err.str()); 
        }
        // else column exists and type matches, proceed!
    }

    // Queue data for a table insert operation into the SQLitePrinter internal buffer
    void SQLitePrinter::insert_data(const unsigned int mpirank, const unsigned long pointID, const std::string& col_name, const std::string& col_type, const std::string& data)
    {
        // Get the pairID for this rank/pointID combination
        std::size_t rowID = pairfunc(mpirank,pointID);
       
        // Make sure we have a record of this column existing in the output table
        // Create it if needed.  
        ensure_column_exists(col_name, col_type);

        // Check if this column exists in the current output buffer
        // Create it if needed
        auto it=buffer_info.find(col_name);
        if(it==buffer_info.end())
        {
            // Column doesn't exist in buffer. Add it.
            std::size_t next_col_index = buffer_info.size();
            buffer_info[col_name] = std::make_pair(next_col_index,col_type);
       
            // Add header data
            buffer_header.push_back(col_name);
            if(buffer_info.size()!=buffer_header.size())
            {
                std::stringstream err;
                err<<"Size of buffer_header ("<<buffer_header.size()<<") does not match buffer_info ("<<buffer_info.size()<<"). This is a bug, please report it.";
                printer_error().raise(LOCAL_INFO.err.str());    
            }

            // Add buffer space
            for(auto jt=transaction_data_buffer.begin(); 
                     jt!=transaction_data_buffer.end(); ++jt)
            {
               // Add new empty column to every row
               // Values are null until we add them 
               jt->push_back("null");

               // Make sure size is correct
               if(jt->size()!=buffer_header.size())
               {
                   std::stringstream err;
                   err<<"Size of a row in the transaction_data_buffer ("<<jt->size()<<") does not match buffer_header ("<<buffer_info.size()<<"). This is a bug, please report it.";
                   printer_error().raise(LOCAL_INFO.err.str());    
               }
            }

            // Now point the map iterator to the right place
            it=buffer_info.find(col_name);
        }
        else
        {
            // Column exists in buffer, but we should also make sure the 
            // type is consistent with the new data we are adding
            std::string buffer_col_type = it->second->second;
            if(buffer_col_type != col_type)
            {
                std::stringstream err;
                err<<"Attempted to add data for column '"<<col_name<<"' to SQLitePrinter transaction buffer, but the type of the new data ("<<col_type<<") does not match the type already recorded for this column in the buffer ("<<buffer_col_type<<").";
                printer_error().raise(LOCAL_INFO.err.str());   
            }
        }

        // Check if a row for this data exists in the transaction buffer
        auto buf_it=transaction_data_buffer.find(rowID);
        if(buf_it=transaction_data_buffer.end())
        {
            // Nope, no row yet for this rowID. Add it.
            // Data is set to 'null' until we add some.
            std::size_t current_row_size=buffer_info.size();
            transaction_data_buffer.emplace(rowID,std::vector<std::string>(current_row_size,"null"));
            auto buf_it=transaction_data_buffer.find(rowID);
        }

        // Add the data to the transaction buffer
        std::size_t col_index = it->second->first;
        transaction_data_buffer.at(rowID).at(col_index) = data;

        // If the buffer is full, execute a transaction to write
        // data to disk, and clear the buffer
        if(transaction_data_buffer.size()>=max_buffer_length)
        {
            dump_buffer();         
        }
    }

    // Delete all buffer data and reset all buffer variables
    void SQLitePrinter::reset_buffer()
    {
        buffer_info.clear();
        buffer_header.clear();
        transaction_data_buffer.clear();
    }

    // Execute an SQLite transaction to write the buffer to the output table
    void SQLitePrinter::dump_buffer()
    {
        require_output_ready(); 
        std::stringstream sql;

        sql<<"INSERT INTO "<<table_name<<" (";
        for(auto col_name_it=buffer_header.begin(); col_name_it!=buffer_header.end(); ++col_name_it)
        {
            sql<<(*col_name_it)<<comma_unless_last(col_name_it,buffer_header);
        }
        sql<<") VALUES (";
        for(auto row_it=transaction_data_buffer.begin(),
                 row_it!=transaction_data_buffer.end(); ++row_it)
        {
            for(auto col_it=row_it->begin(); col_it!=row_it->end(); ++col_it)
            {
                sql<<(*col_it)<<comma_unless_last(col_it,(*row_it));   
            }
            sql<<")"<<comma_unless_last(row_it,transaction_data_buffer);
        }
        sql<<")";

        /* Execute SQL statement */
        int rc;
        char *zErrMsg = 0;
        rc = sqlite3_exec(db, sql.str().c_str(), NULL, NULL, &zErrMsg);
  
        if(rc != SQLITE_OK)
        {
            std::stringstream err;
            err<<"Failed to write transaction buffer for SQLitePrinter to database! SQL error was: "<<zErrMsg;
            sqlite3_free(zErrMsg);
            printer_error().raise(LOCAL_INFO.err.str());        
        }

        // Clear all the buffer data
        reset_buffer(); 
    }
     
    /// @{ PRINT FUNCTIONS
    /// Need to define one of these for every type we want to print!

    /// Templatable print functions
    #define PRINT(TYPE,SQLTYPE) _print(TYPE const& value, const std::string& label, const int vID, const uint rank, const ulong pID) \
       { template_print(value,label,vID,rank,pID,SQLTYPE); }
    void SQLitePrinter::PRINT(int      ,"INTEGER")
    void SQLitePrinter::PRINT(uint     ,"INTEGER")
    void SQLitePrinter::PRINT(long     ,"INTEGER")
    void SQLitePrinter::PRINT(ulong    ,"INTEGER")
    void SQLitePrinter::PRINT(longlong ,"INTEGER")
    void SQLitePrinter::PRINT(ulonglong,"INTEGER")
    void SQLitePrinter::PRINT(float    ,"REAL")
    void SQLitePrinter::PRINT(double   ,"REAL")
    #undef PRINT


  }
}
