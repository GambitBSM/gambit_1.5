//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  HDF5 printer version 2 function defintions
//
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Jan
///
///  *********************************************

#include "gambit/Printers/printers/hdf5printer_v2.hpp"

namespace Gambit
{
  namespace Printers
  {

    /// Retrieve name of the dataset we are supposed to access
    std::string HDF5DataSetBase::myname() const { return _myname; }

    /// Retrieve the dataset ID for the currently open dataset 
    hid_t HDF5DataSetBase::get_dset_id() const
    {
        ensure_dataset_is_open();
        return dset_id;
    }

    /// Retrieve the current size of the dataset on disk
    std::size_t get_dset_length()
    {
        ensure_dataset_is_open();
        return dims[0];
    }

    /// Enforce that the dataset must be open for whatever follows (or else an error is thrown)
    void HDF5DataSetBase::ensure_dataset_is_open() const
    {
        if(not is_open)
        {
            std::ostringstream errmsg;
            errmsg << "Error! Dataset (with name: \""<<myname()<<"\") is not open! Code following this check is not permitted to run. This is a bug in HDF5Printer2, please report it.";
            printer_error().raise(LOCAL_INFO, errmsg.str());     
        }
    }

    /// Open an existing dataset
    hid_t HDF5DataSetBase::open_dataset(hid_t location_id, const std::string& name)
    {
       // Open the dataset
       hid_t out_dset_id = H5Dopen2(location_id, name.c_str(), H5P_DEFAULT);
       if(out_dset_id<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error opening existing dataset (with name: \""<<name<<"\") in HDF5 file. H5Dopen2 failed." << std::endl
                 << "You may have a corrupt hdf5 file from a previous run. Try using -r, or deleting the old output.";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       // Get dataspace of the dataset.
       hid_t dspace_id = H5Dget_space(out_dset_id);
       if(dspace_id<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error opening existing dataset (with name: \""<<name<<"\") in HDF5 file. H5Dget_space failed.";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       // Get the number of dimensions in the dataspace.
       int rank = H5Sget_simple_extent_ndims(dspace_id);
       if(rank<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error opening existing dataset (with name: \""<<name<<"\") in HDF5 file. H5Sget_simple_extent_ndims failed.";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }
       if(rank!=DSETRANK)
       {
          std::ostringstream errmsg;
          errmsg << "Error while accessing existing dataset (with name: \""<<name<<"\") in HDF5 file. Rank of dataset ("<<rank<<") does not match the expected rank ("<<DSETRANK<<").";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       // Get the dimension size of each dimension in the dataspace
       // now that we know ndims matches DSETRANK.
       hsize_t dims_out[DSETRANK];
       int ndims = H5Sget_simple_extent_dims(dspace_id, dims_out, NULL);
       if(ndims<0)
       {
          std::ostringstream errmsg;
          errmsg << "Error while accessing existing dataset (with name: \""<<name<<"\") in HDF5 file. Failed to retrieve dataset extents (H5Sget_simple_extent_dims failed).";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }

       // Update parameters to match dataset contents
       // Compute initial dataspace and chunk dimensions
       dims[0] = dims_out[0]; // Set to match existing data
       maxdims[0] = H5S_UNLIMITED; // No upper limit on number of records allowed in dataset
       chunkdims[0] = CHUNKLENGTH;
       slicedims[0] = 1; // Dimensions of a single record in the data space

       return out_dset_id;
    }

    /// Extend dataset by the specified amount
    void HDF5DataSetBase::extend_dset(const std::size_t extend_by)
    {
       ensure_dataset_is_open();
       std::size_t current_length = dims[0];
       std::size_t newlength = current_length + extend_by;
       dims[0] = newlength; 
       herr_t status = H5Dset_extent(get_dset_id(), dims);
       if(status<0)
       {
          std::ostringstream errmsg;
          errmsg << "Failed to extend dataset (with name: \""<<myname()<<"\") from length "<<current_length<<" to length "<<newlength<<"!";
          printer_error().raise(LOCAL_INFO, errmsg.str());
       }
    }

    /// Obtain memory and dataspace identifiers for writing to a hyperslab in the dataset
    std::pair<hid_t,hid_t> HDF5DataSetBase::select_hyperslab(std::size_t offset, std::size_t length) const
    {
        // Make sure that this chunk lies within the dataset extents
        if(offset + length > dset[0])
        {
           std::ostringstream errmsg;
           errmsg << "Error selecting chunk from dataset (with name: \""<<myname()<<") in HDF5 file. Tried to select a hyperslab which extends beyond the dataset extents:" << std::endl;
           errmsg << "  offset = " << offset << std::endl;
           errmsg << "  offset+length = " << length << std::endl;
           errmsg << "  dset[0] = "<< dset[0] << std::endl;
           printer_error().raise(LOCAL_INFO, errmsg.str());
        }

        // Select a hyperslab.
        hid_t dspace_id = H5Dget_space(get_dset_id());
        if(dspace_id<0) 
        {
           std::ostringstream errmsg;
           errmsg << "Error selecting chunk from dataset (with name: \""<<myname()<<"\") in HDF5 file. H5Dget_space failed." << std::endl;
           printer_error().raise(LOCAL_INFO, errmsg.str());
        }

        hsize_t offsets[DSETRANK];
        offsets[0] = offset;

        hsize_t selection_dims[DSETRANK]; // Set same as output chunks, but may have a different length
        selection_dims[0] = length; // Adjust chunk length to input specification

        herr_t err_hs = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offsets, NULL, selection_dims, NULL);        

        if(err_hs<0) 
        {
           std::ostringstream errmsg;
           errmsg << "Error selecting chunk from dataset (with name: \""<<myname()<<"\", offset="<<offset<<", length="<<selection_dims[0]<<") in HDF5 file. H5Sselect_hyperslab failed." << std::endl;
           printer_error().raise(LOCAL_INFO, errmsg.str());
        }

        // Define memory space
        //H5::DataSpace memspace( DSETRANK, this->get_chunkdims() );
        hid_t memspace_id = H5Screate_simple(DSETRANK, selection_dims, NULL);         

        #ifdef HDF5_DEBUG 
        std::cout << "Debug variables:" << std::endl
                  << "  dset[0]           = " << dset[0] << std::endl
                  << "  offsets[0]        = " << offsets[0] << std::endl
                  << "  CHUNKLENGTH       = " << CHUNKLENGTH << std::endl
                  << "  selection_dims[0] = " << selection_dims[0] << std::endl;
        #endif

        return std::make_pair(memspace_id, dspace_id); // Be sure to close these identifiers after using them!
    }


  }
}

