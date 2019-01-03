//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  HDF5 printer version 2
///
///  This version of the HDF5 printer is a near
///  complete rewrite of the HDF5 printer, based
///  on a 'transaction' concept. Writes are 
///  performed on a single output file in
///  'transactions', during which the file is
///  locked and cannot be accessed by other
///  processes. The idea is based on how the 
///  SQLitePrinter works, though much more manual
///  work is required since HDF5 does not natively
///  work via transactions, unlike SQLite.
///
///  P.S. Also now supports dynamic changing of the
///  buffer size!
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Jan
///
///  *********************************************


#ifndef __hdf5printer_v2_hpp__
#define __hdf5printer_v2_hpp__

#include <algorithm>
#include <set>
#include <iterator>

namespace Gambit
{
  namespace Printers
  {

    /// Length of chunks in chunked HDF5 dataset. Affects write/retrieval performance for blocks of data of various sizes.
    /// It is set to an "intermediate" sort of size since that seems to work well enough. 
    static const std::size_t HDF5_CHUNKLENGTH = 100; 

    template<class T>
    std::set set_diff(const std::set<T>& set1, const std::set<T>& set2)
    {
        std::set<int> result;
        std::set_difference(set1.begin(), set1.end(), set2.begin(), set2.end(),
            std::inserter(result, result.end()));
        return result;
    }

    /// Base class for interfacing to a HDF5 dataset
    class HDF5DataSetBase
    {
      private:

         // Dataset and chunk dimension specification arrays
         // We are only using 1D output datasets for simplicity.
         // Values are only valid if 'is_open==true'
         hsize_t  dims     [DSETRANK];
         hsize_t  maxdims  [DSETRANK];
         hsize_t  chunkdims[DSETRANK];
         hsize_t  slicedims[DSETRANK];
         // Note, dims[0] is current size of dataset, so next unused index is equal to dims[0] 

         /// Name of the dataset in the hdf5 file
         std::string _myname;

         /// Flag to let us known if the dataset is open
         bool is_open;

         /// HDF5 dataset identifer
         hid_t dset_id;

      protected:
 
         /// Dimension of output dataset. We are only using 1D datasets for simplicity.
         static const std::size_t DSETRANK = 1; 
 
         /// Retrieve name of the dataset we are supposed to access
         std::string myname() const;

         /// Enforce that the dataset must be open for whatever follows (or else an error is thrown)
         void ensure_dataset_is_open() const;

         /// Retrieve the dataset ID for the currently open dataset 
         hid_t get_dset_id() const;

         /// Retrieve the current size of the dataset on disk
         std::size_t get_dset_length;

         /// Open dataset on disk and obtain HDF5 handles
         void open_dataset();

         /// Close dataset on disk and release handles
         void close_dataset();

         /// Extend dataset by the specified amount
         void extend_dset(const std::size_t extend_by);
 
         /// Obtain memory and dataspace identifiers for writing to a hyperslab in the dataset
         std::pair<hid_t,hid_t> select_hyperslab(std::size_t offset, std::size_t length) const;
 
    }
 
    /// Class for interfacing to a HDF5 dataset of fixed type
    template<class T>
    class HDF5DataSet: public HDF5DataSetBase
    {
      public:

         /// Write a vector of data to disk at the end of the dataset
         std::size_t write_vector(const std::vector<T>& data)
         {
             open_dataset();
             bool all_data_written=false;
             T buffer[MAX_BUFFER_SIZE];
             std::size_t i = 0;
             std::size_t new_dset_size;
             while(not all_data_written)
             {
                 // Copy data into buffer up to MAX_BUFFER_SIZE
                 for(std::size_t j=0; 
                     (j<MAX_BUFFER_SIZE) && (i<data.size());
                     ++j, ++i)
                 {
                     buffer[j] = data[i];
                 }
                 // Write buffer to disk
                 write_buffer(buffer,j+1);
                 if(i+1==data.size()) all_data_written = true;
             }
             new_dset_size = get_dset_length();
             close_dataset();

             // Report new size of the dataset so that we can check that all datasets are the same length
             return new_dset_size;
         }

         /// Write a block of data to disk at the end of the dataset
         /// This is the lower-level function. There is a fixed-size
         /// buffer that cannot be exceeded. If more data than
         /// MAX_BUFFER_SIZE is to be written then the 'write_vector'
         /// function will split it up and write it in pieces.
         void write_buffer(const T (&buffer)[MAX_BUFFER_SIZE], const std::size_t length) 
         {
             if(length>MAX_BUFFER_SIZE)
             {
                 std::ostringstream errmsg;
                 errmsg << "Error! Received buffer with length ("<<length<<") greater than MAX_BUFFER_SIZE ("<<MAX_BUFFER_SIZE<<") while tring to perform block write for dataset (name="<<myname()<<"). The input to this function is therefore invalid."; 
                 printer_error().raise(LOCAL_INFO, errmsg.str());
             }
             if(length==0)
             {
                 std::ostringstream errmsg;
                 errmsg << "Error! Received buffer of length zero! This will cause an error when trying to select element for writing, and there is no point calling this function with no points to write anyway. Please review the input to this function (error occurred while tring to perform block write for dataset (name="<<myname()<<"))"; 
                 printer_error().raise(LOCAL_INFO, errmsg.str());
             }

             ensure_dataset_is_open();

             // Get the C interface identifier for the type of the output dataset
             hid_t expected_dtype = hdftype_id;
             hid_t dtype = H5Dget_type(get_dset_id()); // type with which the dset was created
             if(not H5Tequal(dtype, expected_dtype))
             {
                 std::ostringstream errmsg;
                 errmsg << "Error! Tried to write to dataset (name="<<myname()<<") with type id "<<dtype<<" but expected it to have type id "<<expected_dtype<<". This is a bug, please report it."; 
                 printer_error().raise(LOCAL_INFO, errmsg.str());
             }

             // Get next unused dataset index
             std::size_t offset = get_dset_length();

             // Extend dataset by the size of the new data to be added
             extend_dset(length);

             // Select output hyperslab
             // (this also determines what data will be read out of the buffer)
             std::pair<hid_t,hid_t> selection_ids = select_hyperslab(offset,length);
             hid_t memspace_id = selection_ids.first;
             hid_t dspace_id   = selection_ids.second;

             // Write the data to the hyperslab.
             herr_t status = H5Dwrite(get_dset_id(), hdftype_id, memspace_id, dspace_id, H5P_DEFAULT, buffer);
             if(status<0)
             {
                std::ostringstream errmsg;
                errmsg << "Error writing new chunk to dataset (with name: \""<<myname()<<"\") in HDF5 file. H5Dwrite failed." << std::endl;
                printer_error().raise(LOCAL_INFO, errmsg.str());
             }
             
             // Release the hyperslab IDs
             H5Sclose(dspace_id);
             H5Sclose(memspace_id);
         }

         /// Write data to disk at specified positions
         void write_random(const std::map<std::size_t,T>& data)
         {
             open_dataset();

             close_dataset();
         }

         /// Write a buffer of data to disk at the specified positions (must be within current dataset extents)
         void write_RA_buffer(const T (&buffer)[MAX_BUFFER_SIZE], const hsize_t (&coords)[MAX_BUFFER_SIZE], std::size_t npoints) 
         {
            if(npoints>MAX_BUFFER_SIZE)
            {
                std::ostringstream errmsg;
                errmsg << "Error! Received npoints ("<<npoints<<") greater than MAX_BUFFER_SIZE ("<<MAX_BUFFER_SIZE<<") while tring to perform RA write for dataset (name="<<myname()<<"). The input to this function is therefore invalid."; 
                printer_error().raise(LOCAL_INFO, errmsg.str());
            }
            if(npoints==0)
            {
                std::ostringstream errmsg;
                errmsg << "Error! Received npoints=0! This will cause an error when trying to select element for writing, and there is no point calling this function with no points to write anyway. Please review the input to this function (error occurred while tring to perform RA write for dataset (name="<<myname()<<"))"; 
                printer_error().raise(LOCAL_INFO, errmsg.str());
            }

            ensure_dataset_is_open();

            bool error_occurred = false; // simple error flag

            // Check that no data is to be written outside the current dataset extents. This
            // function is only for writing back to points that already exist!  
            std::size_t max_coord = *std::max_element(coords,coords+npoints);
            if(max_coord > get_dset_length())
            {
                std::ostringstream errmsg;
                errmsg<<"Attempted to perform RA write to a point outside the current dataset extents! The dataset should be resized prior to calling this function, so this is a bug, please report it.";
                printer_error().raise(LOCAL_INFO, errmsg.str()); 
            }

            // Dataset size in memory
            static const std::size_t MDIM_RANK = 1; 
            hsize_t mdim[] = {npoints};
            
            // Dataspace for the output values
            hid_t dspace = H5Screate_simple(MDIM_RANK, mdim, NULL);
            if(dspace<0) error_occurred = true; 

            // Get the C interface identifier for a copy of the dataspace
            // of the dataset
            hid_t dspace_id = H5Dget_space(get_dset_id());
            if(dspace_id<0) error_occurred = true; 

            // Select the target write points in the file dataspace
            hid_t errflag = H5Sselect_elements(dspace_id, H5S_SELECT_SET, npoints, coords);
            if(errflag<0) error_occurred = true; 

            // Get the C interface identifier for the type of the output dataset
            hid_t expected_dtype = hdftype_id;
            hid_t dtype = H5Dget_type(get_dset_id()); // type with which the dset was created
            if(not H5Tequal(dtype, expected_dtype))
            {
                std::ostringstream errmsg;
                errmsg << "Error! Tried to write to dataset (name="<<myname()<<") with type id "<<dtype<<" but expected it to have type id "<<expected_dtype<<". This is a bug, please report it."; 
                printer_error().raise(LOCAL_INFO, errmsg.str());
            }

            // Write data to selected points
            // (H5P_DEFAULT specifies some transfer properties for the I/O 
            //  operation. These are the default values, probably are ok.)
            hid_t errflag2 = H5Dwrite(get_dset_id(), dtype, dspace, dspace_id, H5P_DEFAULT, values);

            if(errflag2<0) error_occurred = true; 
 
            if(error_occurred)
            {
                std::ostringstream errmsg;
                errmsg << "Error! Failed to write desynchronised buffer data to file! (dataset name="<<myname()<<")"<<std::endl
                       << "Error flags were:" << std::endl
                       << "  dspace   : " << dspace << std::endl
                       << "  dspace_id: " << dspace_id << std::endl
                       << "  errflag  : " << errflag << std::endl
                       << "  errflag2 : " << errflag2 << std::endl
                       << "Variables:" << std::endl
                       << "  dtype = " << dtype;
                printer_error().raise(LOCAL_INFO, errmsg.str());
            }
            
            H5Tclose(dtype);
            H5Sclose(dspace_id);
            H5Sclose(dspace);
         }


      private:

         /// Maximum size of data to be written to dataset in one operation
         static const std::size_t MAX_BUFFER_SIZE = 10000;

         /// HDF5 identifier for the template type of this dataset
         const hid_t hdftype_id;

         /// Create a new dataset at the specified location
         void create_dataset(hid_t location_id, const std::string& name);

    }

    /// HDF5 identifier for the template type of this dataset
    template<class T>
    const hid_t HDF5DataSet::hdftype_id = get_hdf5_data_type<T>::type();
 
    /// Create a (chunked) dataset
    template<class T>
    hid_t HDF5DataSet<T>::create_dataset(hid_t location_id)
    {
       hsize_t  dims     [DSETRANK];
       hsize_t  maxdims  [DSETRANK];
       hsize_t  chunkdims[DSETRANK];
       hsize_t  slicedims[DSETRANK]; 

       // Compute initial dataspace and chunk dimensions
       dims[0] = 0; // Empty to start
       maxdims[0] = H5S_UNLIMITED; // No upper limit on number of records allowed in dataset
       chunkdims[0] = HDF5_CHUNKLENGTH;
       slicedims[0] = 1; // Dimensions of a single record in the data space
 
       // Create the data space
       hid_t dspace_id = H5Screate_simple(DSETRANK, dims, maxdims);
       if(dspace_id<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error creating dataset (with name: \""<<myname()<<"\") in HDF5 file. H5Screate_simple failed.";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       // Object containing dataset creation parameters
       hid_t cparms_id = H5Pcreate(H5P_DATASET_CREATE);
       if(cparms_id<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error creating dataset (with name: \""<<myname()<<"\") in HDF5 file. H5Pcreate failed.";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       herr_t status = H5Pset_chunk(cparms_id, DSETRANK, chunkdims);
       if(status<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error creating dataset (with name: \""<<myname()<<"\") in HDF5 file. H5Pset_chunk failed.";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       // Check if location id is invalid
       if(location_id==-1)
       {
          std::ostringstream errmsg;
          errmsg << "Error! Tried to create hdf5 dataset (with name: \""<<myname()<<"\") at undefined location (location_id was -1). Please check that calling code supplied a valid location handle. This is a bug, please report it.";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       // Create the dataset
       hid_t output_dset_id;
       output_dset_id = H5Dcreate2(location_id, name.c_str(), hdftype_id, dspace_id, H5P_DEFAULT, cparms_id, H5P_DEFAULT);
       if(output_dset_id<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error creating dataset (with name: \""<<myname()<<"\") in HDF5 file. Dataset with same name may already exist";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }
       return output_dset_id;
    }


    /// Base class for buffers
    class HDF5BufferBase
    {
      public:

        /// Make sure buffer includes the input point (data will be set as 'invalid' unless given elsewhere)
        virtual void update(const unsigned int mpirank, const unsigned long pointID) = 0;
 
        /// Empty buffer to disk
        virtual void flush() = 0;

        /// Clear all data in memory ***and on disk*** for this buffer
        virtual void reset() = 0;

    }

    /// Class to manage buffer for a single output label
    template<class T>
    class HDF5Buffer: public HDF5BufferBase
    {
      public:

        /// Make sure buffer includes the specified point (data will be set as 'invalid' unless given elsewhere)
        void update(const unsigned int mpirank, const unsigned long pointID)
        {
            PPIDpair ppid(mpirank,pointID);
            buffer[ppid]; // Create point with default value if it doesn't exist
            buffer_set.insert(ppid);
            auto it = buffer_valid.find(ppid);
            // If point not already in the buffer, set it as invalid
            if(it==buffer_valid.end()) buffer_valid[ppid] = 0;
        }

        /// Insert data to print buffer at the specified point (overwrite if it already exists in the buffer)
        void append(T const& value, unsigned int mpirank, const unsigned long pointID)
        {
            PPIDpair ppid(mpirank,pointID);
            buffer      [ppid] = value;
            buffer_valid[ppid] = 1;
            buffer_set.insert(ppid);
        }

        /// Empty the buffer to disk as block with the specified order
        void block_flush(const std::vector<PPIDpair>& order)
        {
            // Make sure output order is same size as the buffer to be output
            if(order.size() != buffer.size())
            {
                std::ostringstream errmsg;
                errmsg << "Supplied buffer ordering vector is not the same size as the buffer! This is a bug, please report it.";
                printer_error().raise(LOCAL_INFO, errmsg.str());
            }

            // Need to keep track of whether buffer points have been added to the ordered output
            std::set<PPIDpair> done;

            // Create a vector version of the buffer in the specified order
            std::vector<T> ordered_buffer;
            for(auto ppid_it=order.begin(); ppid_it!=order.end(); ++ppid_it)
            {
                if(done.count(*ppid_it)==1)
                {
                    std::ostringstream errmsg;
                    errmsg << "Supplied buffer ordering vector contains a duplicate PPIDpair! This is a bug, please report it.";
                    printer_error().raise(LOCAL_INFO, errmsg.str());
                }
                ordered_buffer.push_back(*ppid_it);
                done.insert(*ppid_it);
            }

            // Check if any points were not added to the ordered buffer
            std::set<PPIDpair> not_done = set_diff(buffer_set,done);
            
            if(not_done.size()>0)
            {
                std::ostringstream errmsg;
                errmsg << "Supplied buffer ordering vector does not specify order positions for all points in the buffer! This is a bug, please report it.";
                printer_error().raise(LOCAL_INFO, errmsg.str());
            }

            // Perform dataset writes
            my_dataset      .write_vector(ordered_buffer);
            my_dataset_valid.write_vector(ordered_buffer);
         }

        /// Empty the buffer to disk as "random access" data at pre-existing positions matching the point IDs
        /// May not completely empty the buffer; points will be removed from the buffer if they are included
        /// in the supplied position map.
        void random_flush(const std::map<PPIDpair,std::size_t>& position_map)
        {
            const std::map<std::size_t,T> pos_buffer;
            const std::map<std::size_t,int> pos_buffer_valid;
            for(auto it=position_map.begin(); it!=position_map.end(); ++it)
            {
                const PPIDpair& ppid = it->first;
                const std::size_t& position = it->second;
                pos_buffer      [position] = buffer      .at(ppid); // TODO: check for duplicates in position_map  
                pos_buffer_valid[position] = buffer_valid.at(ppid); 
            }
          
            // Perform dataset writes          
            my_dataset      .write_random(pos_buffer      );
            my_dataset_valid.write_random(pos_buffer_valid);
         }

        /// Clear all data in the buffer ***and on disk***
        void reset()
        {

        }
      
      private:

        /// Flag to tell us whether this buffer should perform block writes
        /// to the output dataset, or look up and overwrite existing points.
        bool synchronised;

        /// Object that provides an interface to the output HDF5 dataset matching this buffer
        HDF5DataSet<T> my_dataset;
        HDF5DataSet<T> my_dataset_valid;

        /// Buffer containing points to be written to disk upon "flush"
        std::map<PPIDpair,T> buffer;

        /// Set detailing what points are in the buffer
        std::set<PPIDpair> buffer_set; 

        /// Buffer specifying whether the data in the primary buffer is "valid".
        std::map<PPIDpair,int> buffer_valid;

    }

    /// Class to manage a set of buffers for a single output type
    template<class T>
    class HDF5MasterBufferT
    {
      public:

        HDF5Buffer<T>& get_buffer(const std::string& label)
        {
            auto it=my_buffers.find(label);
            if(it==my_buffers.end())
            {
                // No buffer with this name. Need to create one!
                my_buffers.emplace(label,HDF5Buffer<T>());
                it=my_buffers.find(label);
            }
            return (*it);
        }

      private:
 
        std::map<std::string,HDF5Buffer<T>> my_buffers;

    }

    /// Class to manage all buffers for a given printer object
    /// Also handles the file locking/access to the output file
    class HDF5MasterBuffer
    {

      public:

        /// Queue up data to be written to disk when buffers are full
        template<class T>
        void schedule_print(T const& value, const std::string& label, const unsigned int mpirank, const unsigned long pointID)
        {
            HDF5Buffer<T> buffer = get_buffer(label);
            buffer.append(value,mpirank,pointID);

            /// Check if the point is known to be in the buffers already
            PPIDpair thispoint(mpirank,pointID);
            auto it = buffered_points.find(thispoint);
            if(it==buffered_points.end())
            {
                update_all_buffers(mpirank,pointID);
                buffered_points.push_back(thispoint);
            }
        }

        /// Empty all buffers to disk
        void flush()
        {
            for(auto it=all_buffers.begin(); it!=all_buffers.end(); ++it)
            {
                // We have to tell the buffers what order to write their output.
                it->flush(buffered_points);
            } 
        }

        /// Clear all data in buffers ***and on disk*** for this printer
        void reset()
        {
            for(auto it=all_buffers.begin(); it!=all_buffers.end(); ++it)
            {
                it->reset();
            }
        }

      private:

        /// Map containing pointers to all buffers managed by this object;
        std::map<std::string,HDF5BufferBase*> all_buffers;

        /// Vector of PPIDpairs that are currently stored in the printer buffers
        /// This also defines the order in which points should ultimately be
        /// written to disk (we will tell the buffers what order to print stuff,
        /// they don't track the ordering themselves).
        std::vector<PPIDpair> buffered_points;

        /// Retrieve the buffer for a given output label (and type)
        template<class T>
        HDF5Buffer<T> get_buffer(const std::string& label);

        /// Add base class point for a buffer to master buffer map  
        void update_buffer_map(const std::string& label, HDF5BufferBase& buff)
        {
            auto it=all_buffers.find(label);
            if(it==all_buffers.end())
            {
                /// Not already in the map; add it
                all_buffers.emplace(label,&buff);
            }
        }

        /// Inform all buffers that data has been written to certain mpirank/pointID pair
        /// They will make sure that they have an output slot for this pair, so that all the
        /// buffers for this printer stay "synchronised".
        void update_all_buffers(const unsigned int mpirank, const unsigned long pointID)
        {
            for(auto it=all_buffers.begin(); it!=all_buffers.end(); ++it)
            {
                it->update(mpirank,pointID);
            }
        }

        /// Buffer manager objects
        //  Need a map for every directly printable type, and a specialisation
        //  of 'get_buffer' to access it. But the latter have to be
        //  defined outside the class declaration, so they can be found below.
        HDF5MasterBufferT<int      > hdf5_buffers_int;
        HDF5MasterBufferT<uint     > hdf5_buffers_uint;
        HDF5MasterBufferT<long     > hdf5_buffers_long;
        HDF5MasterBufferT<ulong    > hdf5_buffers_ulong;
        HDF5MasterBufferT<longlong > hdf5_buffers_longlong;
        HDF5MasterBufferT<ulonglong> hdf5_buffers_ulonglong;
        HDF5MasterBufferT<float    > hdf5_buffers_float;
        HDF5MasterBufferT<double   > hdf5_buffers_double;
 
    }

    /// Specialisations for 'get_buffer' function for each buffer type
    #define SPECIALISE_GET_BUFFER(TYPE)\
    template<>\
    HDF5Buffer<TYPE>& HDF5MasterBuffer::get_buffer<TYPE>(const std::string& label)\
    {\
        HDF5Buffer<TYPE>& out_buffer = CAT(hdf5_buffers_,TYPE).get_buffer(label);\
        update_buffer_map(label,out_buffer);\
        return out_buffer;\
    }
    SPECIALISE_GET_BUFFER(int      )
    SPECIALISE_GET_BUFFER(uint     )
    SPECIALISE_GET_BUFFER(long     )
    SPECIALISE_GET_BUFFER(ulong    )
    SPECIALISE_GET_BUFFER(longlong )
    SPECIALISE_GET_BUFFER(ulonglong)
    SPECIALISE_GET_BUFFER(float    )
    SPECIALISE_GET_BUFFER(double   )










    /// The main printer class for output to HDF5 format
    class HDF5Printer2 : public BasePrinter
    {
      public:
        /// Constructor (for construction via inifile options)
        HDF5Printer2(const Options&, BasePrinter* const primary = NULL);

        /// Destructor
        ~HDF5Printer2();

        /// Base class virtual function overloads
        /// (the public virtual interface)
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

      private:

    };




  }
}

#endif
