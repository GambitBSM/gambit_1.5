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
///
///  TODO: Turns out SQLite is case-insensitive, so need
///  to change various string comparisons here to also be
///  case-insensitive.

#include <iostream>
#include <sstream>
#include <chrono>
#include <thread>

// SQLite3 C interface 
#include <sqlite3.h> 

// Gambit
#include "gambit/Printers/printers/sqliteprinter.hpp"

// Define this macro to dump attempted SQL statements during exceptions
#define SQL_DEBUG

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

    // Returns new iterator pointing to next element
    template <typename Iter>
    Iter next_el(Iter iter)
    {
        return ++iter;
    }

    // Return a comma unless iterator is the last in the supplied iterable
    // (or if it points to end())
    template <typename Iter, typename Cont>
    std::string comma_unless_last(Iter it, const Cont& c)
    { 
       std::string out("");
       if((it == c.end()) || (next_el(it) == c.end()))
       { /* this is the last element or end(), do nothing */ }
       else
       { out = ","; }
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
    static int col_name_callback(void* colmap_in, int /*count*/, char** data, char** /* columns */)
    {
        typedef std::map<std::string, std::string, Utils::ci_less> mymaptype;
        mymaptype *colmap = static_cast<mymaptype*>(colmap_in);
   
        // We know that the column name is the second column of the results set, and the
        // type is the third column
        std::string column_name(data[1]);
        std::string column_type(data[2]);

        //std::cout<<"Reading existing columns: name: "<<column_name<<", type: "<<column_type<<std::endl;

        // Add to map
        (*colmap)[column_name] = column_type;
    
        return 0;
    }

    // Constructor
    SQLitePrinter::SQLitePrinter(const Options& options, BasePrinter* const primary)
    : BasePrinter(primary,options.getValueOrDef<bool>(false,"auxilliary"))
#ifdef WITH_MPI
    , myComm() // initially attaches to MPI_COMM_WORLD
#endif
    , mpiRank(0)
    , mpiSize(1)
    , primary_printer(NULL)
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
    , synchronised(!options.getValueOrDef<bool>(false,"auxilliary"))
    {
        // Test options?
        //std::cout << "options:" << std::endl<<options.getNode()<<std::endl;

        if(is_auxilliary_printer())
        {
            // If this is an "auxilliary" printer then we need to get some
            // of our options from the primary printer
            primary_printer = dynamic_cast<SQLitePrinter*>(this->get_primary_printer());
            database_file     = primary_printer->get_database_file(); 
            table_name        = primary_printer->get_table_name();
            max_buffer_length = primary_printer->get_max_buffer_length();
        }
        else
        {
            // MPI setup
#ifdef WITH_MPI 
            this->setRank(myComm.Get_rank()); // tells base class about rank
            mpiRank = myComm.Get_rank();
            mpiSize = myComm.Get_size();
#endif

            // Tell scannerbit if we are resuming
            set_resume(options.getValue<bool>("resume"));

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
            table_name = options.getValueOrDef<std::string>("results","table_name");

            // Delete final target file if one with same name already exists? (and if we are restarting the run)
            // Mostly for convenience during testing. Recommend to use 'false' for serious runs to avoid
            // accidentally deleting valuable output.
            bool overwrite_file  = options.getValueOrDef<bool>(false,"delete_file_on_restart");

            if(getRank()==0 and overwrite_file and not get_resume())
            {
                // Note: "not resume" means "start or restart"
                // Delete existing output file
                std::ostringstream command;
                command << "rm -f "<<database_file;
                logger() << LogTags::printers << LogTags::info << "Running shell command: " << command.str() << EOM;
                FILE* fp = popen(command.str().c_str(), "r");
                if(fp==NULL)
                {
                    // Error running popen
                    std::ostringstream errmsg;
                    errmsg << "rank "<<getRank()<<": Error deleting existing output file (requested by 'delete_file_on_restart' printer option; target filename is "<<database_file<<")! popen failed to run the command (command was '"<<command.str()<<"')";
                    printer_error().raise(LOCAL_INFO, errmsg.str());
                }
                else if(pclose(fp)!=0)
                {
                    // Command returned exit code!=0, or pclose failed
                    std::ostringstream errmsg;
                    errmsg << "rank "<<getRank()<<": Error deleting existing output file (requested by 'delete_file_on_restart' printer option; target filename is "<<database_file<<")! Shell command failed to executed successfully, see stderr (command was '"<<command.str()<<"').";
                    printer_error().raise(LOCAL_INFO, errmsg.str());
                }
            }

#ifdef WITH_MPI
            // Make sure no processes try to open database until we are sure it won't be deleted and replaced
            myComm.Barrier();
#endif
        }       
 
        // Create/open the database file 
        open_db(database_file);

        // Create the results table in the database (if it doesn't already exist)
        make_table(table_name);

        // If we are rank 0 and also resuming (and this is the primary printer), need to read the database and find the previous
        // highest pointID numbers used.
        if(not is_auxilliary_printer())
        { 
            // Record of the highest pointIDs in previous run that are relevant for *this* run
            std::vector<std::size_t> highests(mpiSize);

            if(getRank()==0 and get_resume())
            {
                // Map from ranks to highest pointIDs in previous output
                std::map<std::size_t, std::size_t> highest_pointIDs;
  
                // Construct the SQLite3 statement to retrieve previous ranks and pointIDs from the database
                std::stringstream sql;
                sql << "SELECT MPIrank,pointID FROM "<<table_name;
    
                /* Execute SQL statement and iterate through results*/ 
                sqlite3_stmt *stmt;
                int rc = sqlite3_prepare_v2(db, sql.str().c_str(), -1, &stmt, NULL);
                if (rc != SQLITE_OK) {
                    std::stringstream err;
                    err<<"Encountered SQLite error while preparing to retrieve previous pointIDs: "<<sqlite3_errmsg(db);
                    printer_error().raise(LOCAL_INFO, err.str());
                }
                while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
                    std::size_t rank = sqlite3_column_int(stmt, 0);
                    std::size_t pID  = sqlite3_column_int(stmt, 1);
                    if( pID > highest_pointIDs[rank] )
                    {
                        highest_pointIDs[rank] = pID;
                    }
                }
                if (rc != SQLITE_DONE) {
                    std::stringstream err;
                    err<<"Encountered SQLite error while retrieving previous pointIDs: "<<sqlite3_errmsg(db);
                    printer_error().raise(LOCAL_INFO, err.str());
                }
                sqlite3_finalize(stmt);
 
                // Grab only the highest pointIDs for ranks that we are using THIS run.              
                for (size_t rank = 0; rank < myComm.Get_size(); rank++ )
                {
                    auto it = highest_pointIDs.find(rank);
                    if (it != highest_pointIDs.end())
                        highests[rank] = it->second;
                    else
                        highests[rank] = 0;
                }
            } 

            // Need to communicate these to ScannerBit so it doesn't start assigning them
            // again from zero.
#ifdef WITH_MPI
            int resume_int = get_resume();
            myComm.Barrier();
            // Make sure everyone agrees on the resume status. Not really needed, but I
            // just copied it from the HDF5 printer.
            myComm.Bcast(resume_int, 1, 0);
            set_resume(resume_int);

            if (get_resume())
            {
                std::size_t my_highest;
                myComm.Barrier();
                myComm.Scatter(highests, my_highest, 0);
                get_point_id() = my_highest;
            }
#else
            if (get_resume())
            {
                get_point_id() = highests[0];
            }
#endif
        }
    }

    std::string SQLitePrinter::get_database_file() {return database_file;}
    std::string SQLitePrinter::get_table_name() {return table_name;}
    std::size_t SQLitePrinter::get_max_buffer_length() {return max_buffer_length;}

    void SQLitePrinter::initialise(const std::vector<int>&)
    {
        // Don't need to initialise anything for this printer
    }

    void SQLitePrinter::reset(bool /*force*/)
    {
        // This is needed by e.g. MultiNest to delete old weights and replace them
        // with new ones.

        // Should put a check so that only auxilliary printers can delete stuff, unless 'force' is set

        // Read through header to see what columns this printer has been touching. These are
        // the ones that we will reset/delete.
        // (a more nuanced reset might be required in the future?)
        sql<<"INSERT INTO "<<table<<" (pairID,";
        for(auto col_name_it=buffer_header.begin(); col_name_it!=buffer_header.end(); ++col_name_it)
        {
            sql<<"`"<<(*col_name_it)<<"`"<<comma_unless_last(col_name_it,buffer_header);
        }
        sql<<") VALUES ";
 
    }

    void SQLitePrinter::finalise(bool /*abnormal*/)
    {
        // Dump buffer to disk. Nothing special needed for early shutdown.
        SQLitePrinter::dump_buffer();
    }

    // Reader construction options for constructing a reader
    // object that can read the output we are printing
    Options SQLitePrinter::resume_reader_options()
    {
        Options options;
        // Set options that we need later to construct a reader object for
        // previous output, if required.
        options.setValue("type", "sqlite");
        options.setValue("file", database_file);
        options.setValue("table", table_name);
        return options;
    }

    // Open database and 'attach' it to this object
    void SQLitePrinter::open_db(const std::string& path)
    { 
        // Check if we already have an open database
        if(db!=NULL)
        {
            std::stringstream err;
            err << "Refused to open database file '"<<path<<"'; a database file pointer has already been attached to the SQLite printer!";
            printer_error().raise(LOCAL_INFO,err.str()); 
        }

        if(db_is_open)
        {
            std::stringstream err;
            err << "Refused to open database file '"<<path<<"'; a database file is already flagged by the SQLite printer as open!";
            printer_error().raise(LOCAL_INFO,err.str());  
        }

        int rc; // return code

        rc = sqlite3_open(path.c_str(), &db);
   
        if( rc ) 
        {
            std::stringstream err;
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
            std::stringstream err;
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

    // Function to repeatedly attempt an SQLite statement if the database is locked/busy
    int SQLitePrinter::submit_sql(const std::string& local_info, const std::string& sqlstr, bool allow_fail, sql_callback_fptr callback, void* data, char **zErrMsg)
    {
        int rc;
        do
        {
            rc = sqlite3_exec(db, sqlstr.c_str(), callback, data, zErrMsg);
            if(rc==SQLITE_BUSY)
            {
                // Wait at least a short time to avoid slamming the filesystem too much
                std::chrono::milliseconds timespan(10);
                std::this_thread::sleep_for(timespan);
            }
        }
        while (rc == SQLITE_BUSY);

        // if allow_fail is true then we don't catch this error, we allow the caller of this function to hander it.
        if( (rc != SQLITE_OK) and not allow_fail ){
            std::stringstream err;
            err << "SQL error: " << *zErrMsg << std::endl;
#ifdef SQL_DEBUG
            err << "The attempted SQL statement was:"<<std::endl;
            err << sqlstr << std::endl;; 
#endif
            sqlite3_free(*zErrMsg);
            printer_error().raise(local_info,err.str());
       }
       return rc;
    }

    // Create results table
    void SQLitePrinter::make_table(const std::string& name)
    {
        // Construct the SQLite3 statement
        std::stringstream sql;
        sql << "CREATE TABLE IF NOT EXISTS "<<name<<"("
            << "pairID   INT PRIMARY KEY NOT NULL,"
            << "MPIrank  INT,"
            << "pointID  INT"
            << ");";

        /* Execute SQL statement */
        submit_sql(LOCAL_INFO, sql.str());

        // Flag the results table as existing
        results_table_exists = true;
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
            sql<<"ALTER TABLE "<<table_name<<" ADD COLUMN `"<<sql_col_name<<"` "<<sql_col_type<<";";

            /* Execute SQL statement */
            int rc;
            char *zErrMsg = 0;
            // Need allow_fail=true for this case
            rc = submit_sql(LOCAL_INFO, sql.str(), true, NULL, NULL, &zErrMsg);
  
            if( rc != SQLITE_OK ){
                // Operation failed for some reason. Probably because the column already
                // exists, but we better make sure.

                std::stringstream sql2;
                sql2<<"PRAGMA table_info("<<table_name<<");";

                /* Execute SQL statement */ 
                int rc2;
                char *zErrMsg2 = 0;
                std::map<std::string, std::string, Utils::ci_less> colnames; // Will be passed to and filled by the callback function
                rc2 = submit_sql(LOCAL_INFO, sql2.str(), true, &col_name_callback, &colnames, &zErrMsg2);
 
                if( rc2 != SQLITE_OK ){
                    std::stringstream err;
                    err << "Failed to check SQL column names in output table, after failing to add a new column '"<<sql_col_name<<"' to that table."<<std::endl; 
                    err << "  First SQL error was: " << zErrMsg << std::endl;
#ifdef SQL_DEBUG
                    err << "  The attempted SQL statement was:"<<std::endl;
                    err << sql.str() << std::endl; 
#endif
                    err << "  Second SQL error was: " << zErrMsg2 << std::endl;
#ifdef SQL_DEBUG
                    err << "  The attempted SQL statement was:"<<std::endl;
                    err << sql2.str() << std::endl; 
#endif
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
#ifdef SQL_DEBUG
                    err << "The attempted SQL statement was:"<<std::endl;
                    err << sql.str() << std::endl;
#endif
                    sqlite3_free(zErrMsg);
                    printer_error().raise(LOCAL_INFO,err.str());
                }
                else if(!Utils::iequals(jt->second,sql_col_type))
                {
                    // NOTE: All sorts of type names are equivalent, so this simple string checking is
                    // totally unreliable! 

                    // // Column found, but has the wrong type
                    // std::stringstream err;
                    // err << "Failed to add new column '"<<sql_col_name<<"' to output SQL table! The column already exists, but it has the wrong type (existing column has type '"<<jt->second<<"', but we expected it to have type '"<<sql_col_type<<"'!";
                    // sqlite3_free(zErrMsg);
                    // printer_error().raise(LOCAL_INFO,err.str());
                }

                // Column exists and has the right type! So everything is ok after all.
            }

            // Column should exist now. Need to add the fact of this columns existence to our internal record.
            column_record[sql_col_name] = sql_col_type;     
        } 
        else if(!Utils::iequals(it->second,sql_col_type))
        {
            // // Records say column exists, but not with the type requested!
            // NOTE: All sorts of type names are equivalent, so this simple string checking is
            // totally unreliable! 

            // std::stringstream err;
            // err << "SQLitePrinter records indicated that the column '"<<sql_col_name<<"' already exists in the output table, but with a different type than has been requested (existing type is '"<<it->second<<"', requested type was '"<<sql_col_type<<"'). This indicates either duplicate names in the printer output, or an inconsistency in how the print commands have been issued.";
            // printer_error().raise(LOCAL_INFO,err.str()); 
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

        // Check if a row for this data exists in the transaction buffer
        auto buf_it=transaction_data_buffer.find(rowID);
        if(buf_it==transaction_data_buffer.end())
        {
            // Nope, no row yet for this rowID. Add it.
            // But we should first dump the buffer if it was full
        
            // If the buffer is full, execute a transaction to write
            // data to disk, and clear the buffer
            if(transaction_data_buffer.size()>=max_buffer_length)
            {
                dump_buffer();         
            }
    
            // Data is set to 'null' until we add some.
            std::size_t current_row_size=buffer_info.size();
            transaction_data_buffer.emplace(rowID,std::vector<std::string>(current_row_size,"null"));
        }

        // Check if this column exists in the current output buffer
        // Create it if needed
        auto it=buffer_info.find(col_name);
        if(it==buffer_info.end())
        {
            // Column doesn't exist in buffer. Add it.
            std::size_t next_col_index = buffer_info.size();
            buffer_info[col_name] = std::make_pair(next_col_index,col_type);
       
            // Add header data
            //std::cout<<"Adding column to buffer: "<<col_name<<std::endl;
            buffer_header.push_back(col_name);
            if(buffer_info.size()!=buffer_header.size())
            {
                std::stringstream err;
                err<<"Size of buffer_header ("<<buffer_header.size()<<") does not match buffer_info ("<<buffer_info.size()<<"). This is a bug, please report it.";
                printer_error().raise(LOCAL_INFO,err.str());    
            }

            // Add buffer space
            for(auto jt=transaction_data_buffer.begin(); 
                     jt!=transaction_data_buffer.end(); ++jt)
            {
               std::vector<std::string>& row = jt->second;

               // Add new empty column to every row
               // Values are null until we add them 
               row.push_back("null");

               // Make sure size is correct
               if(row.size()!=buffer_header.size())
               {
                   std::stringstream err;
                   err<<"Size of a row in the transaction_data_buffer ("<<row.size()<<") does not match buffer_header ("<<buffer_info.size()<<"). This is a bug, please report it.";
                   printer_error().raise(LOCAL_INFO,err.str());    
               }
            }

            // Now point the map iterator to the right place

            it=buffer_info.find(col_name);
        }
        else
        {
            // Column exists in buffer, but we should also make sure the 
            // type is consistent with the new data we are adding
            std::string buffer_col_type = it->second.second;
            if(!Utils::iequals(buffer_col_type,col_type))
            {
                std::stringstream err;
                err<<"Attempted to add data for column '"<<col_name<<"' to SQLitePrinter transaction buffer, but the type of the new data ("<<col_type<<") does not match the type already recorded for this column in the buffer ("<<buffer_col_type<<").";
                printer_error().raise(LOCAL_INFO,err.str());   
            }
        }

        // Add the data to the transaction buffer
        std::size_t col_index = it->second.first;
        transaction_data_buffer.at(rowID).at(col_index) = data;
    }

    // Delete all buffer data. Leaves the header intact so that we know what columns
    // this printer has been working with (needed so we can reset them if needed!)
    void SQLitePrinter::clear_buffer()
    {
        transaction_data_buffer.clear();
    }

    // Create an SQL table insert operation for the current transaction_data_buffer
    // Modifies 'sql' stringstream in-place
    void SQLitePrinter::turn_buffer_into_insert(std::stringstream& sql, const std::string& table)
    {
        sql<<"INSERT INTO "<<table<<" (pairID,";
        for(auto col_name_it=buffer_header.begin(); col_name_it!=buffer_header.end(); ++col_name_it)
        {
            sql<<"`"<<(*col_name_it)<<"`"<<comma_unless_last(col_name_it,buffer_header);
        }
        sql<<") VALUES ";
        for(auto row_it=transaction_data_buffer.begin();
                 row_it!=transaction_data_buffer.end(); ++row_it)
        {
            sql<<"(";
            std::size_t pairID = row_it->first;
            std::vector<std::string>& row = row_it->second;
            sql<<pairID<<",";
            for(auto col_it=row.begin(); col_it!=row.end(); ++col_it)
            {
                sql<<(*col_it)<<comma_unless_last(col_it,row);   
            }
            sql<<")"<<comma_unless_last(row_it,transaction_data_buffer);
        }
        sql<<";"; // End statement
    }

    // Execute an SQLite transaction to write the buffer to the output table
    void SQLitePrinter::dump_buffer_as_INSERT()
    {
        // Add the table INSERT operation to a stream
        std::stringstream sql;
        turn_buffer_into_insert(sql,table_name);

        /* Execute SQL statement */
        submit_sql(LOCAL_INFO,sql.str());
    }
 
    void SQLitePrinter::dump_buffer_as_UPDATE()
    {
        std::stringstream sql;
        // So for this is seems like the best thing to do is create a temporary
        // table with this new data, and then update the main output table from
        // this. Otherwise we have to write tonnes of separate 'update' statements,
        // which is probably not very fast.
        // So first we need to create the temporary table.
        sql << "DROP TABLE IF EXISTS temp_table;"
            << "CREATE TEMPORARY TABLE temp_table("
            << "pairID   INT PRIMARY KEY NOT NULL,";
        for(auto col_it=buffer_info.begin(); col_it!=buffer_info.end(); ++col_it)
        { 
            const std::string& col_name(col_it->first);
            const std::string& col_type(col_it->second.second);
            sql<<"`"<<col_name<<"`   "<<col_type<<comma_unless_last(col_it,buffer_info);
        }
        sql <<");";

        // Insert data into the temporary table
        turn_buffer_into_insert(sql,"temp_table");

        // Update the primary output table using the temporary table
        // Following: https://stackoverflow.com/a/47753166/1447953
        sql<<"UPDATE "<<table_name<<" SET (";
        for(auto col_name_it=buffer_header.begin(); col_name_it!=buffer_header.end(); ++col_name_it)
        { 
            sql<<*col_name_it<<comma_unless_last(col_name_it,buffer_header);
        }
        sql<<") = (SELECT ";
        for(auto col_name_it=buffer_header.begin(); col_name_it!=buffer_header.end(); ++col_name_it)
        { 
            sql<<"temp_table."<<*col_name_it<<comma_unless_last(col_name_it,buffer_header);
        }
        sql<<" FROM temp_table WHERE temp_table.pairID = "<<table_name<<".pairID)";
        sql<<" WHERE EXISTS ( SELECT * FROM temp_table WHERE temp_table.pairID = "<<table_name<<".pairID);";

        /* Execute SQL statement */
        submit_sql(LOCAL_INFO,sql.str());
    }
   
    void SQLitePrinter::dump_buffer()
    {
        require_output_ready(); 
        if(synchronised)
        {
            // Primary dataset writes can be performed as INSERT operations
            dump_buffer_as_INSERT();
        }
        else
        {
            // Asynchronous ('auxilliary') writes need to be performed as UPDATE operations
            dump_buffer_as_UPDATE();
        }
        // Clear all the buffer data
        clear_buffer(); 
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

    void SQLitePrinter::_print(ModelParameters const& value, const std::string& label, const int vID, const unsigned int mpirank, const unsigned long pointID)
    {
      std::map<std::string, double> parameter_map = value.getValues();
      _print(parameter_map, label, vID, mpirank, pointID);
    }

    void SQLitePrinter::_print(const map_str_dbl& map, const std::string& label, const int vID, const unsigned int mpirank, const unsigned long pointID)
    {
      for (std::map<std::string, double>::const_iterator
           it = map.begin(); it != map.end(); it++)
      {
        std::stringstream ss;
        ss<<label<<"::"<<it->first;
        _print(it->second, ss.str(), vID, mpirank, pointID);
      }
    }

    void SQLitePrinter::_print(std::vector<double> const& value, const std::string& label, const int vID, const unsigned int mpirank, const unsigned long pointID)
    {
      for(unsigned int i=0;i<value.size();i++)
      {
        std::stringstream ss;
        ss<<label<<"["<<i<<"]";
        _print(value.at(i), ss.str(), vID, mpirank, pointID);  
      }
    }



  }
}
