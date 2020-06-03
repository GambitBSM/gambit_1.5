//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Declaration of the Classy_input class used for
///  communicating with the backend classy.
///  It has the attribute 'input_dict' which is a python
///  dictionary containing the input parameters for CLASS
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2018 May
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///  \date 2019 Mar
///
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///  \date 2020 Mar
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2020 Apr
///
///  *********************************************

#pragma once

#include "gambit/Utils/util_types.hpp"

#ifdef HAVE_PYBIND11

  #include <pybind11/stl.h>

  namespace Gambit
  {

    /// Merge two dictionaries containing input parameters for CLASS
    /// (specific rules are applied in case of key duplication)
    void merge_pybind_dicts(pybind11::dict&, pybind11::dict, bool);

    /// helper to convert the memory address a double pointer points to
    /// to an integer (-> uintptr_t, size of type depends on system & ensures
    /// it is big enough to store memory addresses of the underlying setup)
    /// (implemented to pass contents of arrays to CLASS)
    uintptr_t memaddress_to_uint(double* ptr);

    // Class that manages the input dictionary for classy
    class Classy_input
    {
      public:

        /// add all entries from extra_entries to input_dict, concatenates and returns all
        /// keys that are contained in both dictionaries:
        /// -> no keys in common: returns empty string ("")
        /// -> else: returns sting containing all duplicated keys
        /// need to check after use of this function if returned string was empty to avoid overwriting of
        /// input values & inconsistencies.
        std::string add_dict(pybind11::dict);

        void add_entry(str, double);
        void add_entry(str, int);
        void add_entry(str, str);
        void add_entry(str, std::vector<double>&);

        bool has_key(str);

        // merge dictionaries with overwriting/combining rules that only
        // apply for CLASS input dictionaries
        void merge_input_dicts(pybind11::dict);

        // routine to print CLASS input values to logger
        std::string print_entries_to_logger();

        // clears all entries from input_dict
        void clear();

       // return input_dict
        pybind11::dict get_input_dict();

      private:
        pybind11::dict input_dict;
    };

  }

#endif
