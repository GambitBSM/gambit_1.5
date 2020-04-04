//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of the Classy_input class used for
///  communicating with the backend classy.
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

#include "gambit/Backends/backend_types/classy.hpp"

namespace Gambit
{

  /// helper to convert the memory address a double pointer points to
  /// to an integer (-> uintptr_t, size of type depends on system & ensures
  /// it is big enough to store memory addresses of the underlying setup)
  /// (implemented to pass contents of arrays to CLASS)
  uintptr_t memaddress_to_uint(double* ptr)
  {
    uintptr_t addr;
    addr = reinterpret_cast<uintptr_t>(ptr);
    return addr;
  }

  /// function to merge python dictionary extra_dict into input_dict. If both dictionaries have the same key
  /// it has to be decided on a case-by-case basis how to deal with that
  /// some examples are already implemented, might have to be extended in the future
  /// if more CLASS features are used.
  /// -> very specific for merging the classy input dictionaries
  /// typical example for this would be the key 'output':
  ///    dict input_dict: 'output' : 'tCl nCl' and
  ///    dict extra_dict: 'output' : 'tCl mPk' => mPk has to be added
  ///    => results in    'output' : 'tCl nCl mPk'
  void merge_pybind_dicts(pybind11::dict& input_dict, pybind11::dict extra_dict, bool first_run)
  {

    // loop through extra_dict (should typically be the shorter one)
    for (auto item : extra_dict)
    {
      pybind11::str key = pybind11::str(item.first);
      pybind11::str arg = pybind11::str(item.second);

      // if item not contained in extra_dict but not in input_dict it will be added to input_dict
      if(!input_dict.attr("has_key")(key).cast<bool>())
      {
        input_dict[key] = arg;
        //std::cout << "Adding key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
      }

      // if item contained in both: have to decide on a case-by-case basis what to do!
      else
      {
        // if 'output' is defined twice the entries have to be merged
        // e.g. input_dict['output'] = 'A B C'
        //      extra_dict['output'] = 'A B X Y'
        // should result in input_dict['output'] = 'A B C X Y'
        // (python string.find("x") returns -1 if "x" not contained)
        if(key.attr("find")(pybind11::str("output")).cast<int>()!=-1 ||
           key.attr("find")(pybind11::str("modes")).cast<int>()!=-1)
        {
          // split string of extra_dict['output'] by spaces into list
          // (-> in the example above this would give list = ['A', 'B', 'X', 'Y'])
          pybind11::list list = extra_dict[key].attr("split")();
          // iterate through list and check if current item is also contained in input_dict['output']
          // if it is not it will be added
          for(auto it : list)
          {
            std::string list_entry = std::string(pybind11::str(it));

            // add entry if it is not already in input_dict['output']
            if(input_dict[key].attr("find")(pybind11::str(list_entry)).cast<int>()==-1)
            {
              // add part of extra_dict[key] string that is not contained in input_dict[key] string to input_dict[key]
              std::string new_arg=std::string(pybind11::str(input_dict[key]))+std::string(" ")+std::string(list_entry);
              input_dict[key]= new_arg;
            }
          }
        }

        // if 'l_max...' (scalars or tensors) is defined twice use the maximum requested value
        // e.g. input_dict['l_max_scalars'] = '500'
        //      extra_dict['l_max_scalars'] = '2500'
        // should result in input_dict['l_max_scalars'] = '2500'
        else if(key.attr("find")(pybind11::str("l_max")).cast<int>()!=-1)
        {
          // cast pybind11::detail::item_accessor to pybind11::str, to c++ string and then to int
          //  (I know... the problem is that the yaml file entries are all parsed as strings so
          //  we don't have to distinguish between the different types of different CLASS input keys.
          //  The python wrapper parses everything as string so it is fine. Only here, where we actually
          //  have to compare the passed numbers this causes troubles and we have to cast from a python
          //  string to a c++ int to be able to tell which entry is the larger one. If you have an idea
          //  to make this nicer --> go for it!!)
          int lmax_input_dict = std::stoi(std::string(pybind11::str(input_dict[key])));
          int lmax_extra_dict = std::stoi(std::string(pybind11::str(extra_dict[key])));

          // if lmax_extra_dict is higher than the entry in the input_dict replace it
          if (lmax_input_dict < lmax_extra_dict){input_dict[key] = lmax_extra_dict;}
          // if not 'input_dict'already contains the higher value and there is nothing to do here.
        }
        else
        {
          if(first_run)
          {
            // (JR) TODO see what other combinations could go wrong here -> different likelihoods
            // asking for different redshift bins?
            // But need this check only on the first run (the inputs we are worried about are
            // parameters specifying which output CLASS should produce so they won't change
            // in-between the calculation of different points in one scan)
            std::cout <<"___________________________________________________________________"<<std::endl;
            std::cout <<"___________________________________________________________________"<<std::endl;
            std::cout << "Both dictionaries to merge contain key" << std::string(key) << "with entries "
                << std::string(pybind11::str(input_dict[key]))<< " and " << std::string(pybind11::str(extra_dict[key])) << ". Don't know how to deal with that, yet. Will be taken care of soon." << std::endl;
            std::cout <<"___________________________________________________________________"<<std::endl;
            std::cout <<"___________________________________________________________________"<<std::endl;
          }
        }
      }
    }
  }

  void Classy_input::add_entry(str key, double value) {input_dict[key.c_str()]=std::to_string(value).c_str();}
  void Classy_input::add_entry(str key, int    value) {input_dict[key.c_str()]=std::to_string(value).c_str();}
  void Classy_input::add_entry(str key, str    value) {input_dict[key.c_str()]=value.c_str();}
  void Classy_input::add_entry(str key, std::vector<double>& values)
  {
      // get pointers to arrays holding the information that needs
      // to be passed on to class, convert to uintptr_t (type large enough
      // to store memory address of the used system) and pass to class
      input_dict[key.c_str()] = memaddress_to_uint(&values[0]);
  }

  bool Classy_input::has_key(str key){return input_dict.contains(key.c_str());}

  /// add all entries from extra_entries to input_dict, concatenates and returns all
  /// keys that are contained in both dictionaries:
  /// -> no keys in common: returns empty string ("")
  /// -> else: returns sting containing all duplicated keys
  /// need to check after use of this function if returned string was empty to avoid overwriting of
  /// input values & inconsistencies.
  std::string Classy_input::add_dict(pybind11::dict extra_dict)
  {
    // string to be returned -- stays empty if no duplicated dict entries are found
    std::string common_keys ("");

    for (auto item : extra_dict)
    {
      pybind11::str key = pybind11::str(item.first);
      pybind11::str arg = pybind11::str(item.second);

      if(!input_dict.attr("has_key")(key).cast<bool>())
      {
        input_dict[key] = arg;
        //std::cout << "Adding key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
      }
      else
      {
        // add duplicated strings
        common_keys = common_keys + " " + key.cast<std::string>();
        //std::cout << "Found duplicated key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
      }
    }
    return common_keys;
  }

  // function to merge python dictionary extra_dict into input_dict (member of Classy_input).
  void Classy_input::merge_input_dicts(pybind11::dict extra_dict)
  {
    static bool first_run = true;
    // call function to merge classy pybind dicts defined in CosmoBit_utils.cpp
    // this function takes care of merging the input. If a key is contained in both
    // dictionaries it has to be decided on a case-by-case basis how to deal with that
    // some examples are already implemented, might have to be extended in the future
    // if more CLASS features are used.
    merge_pybind_dicts(input_dict, extra_dict, first_run);
    first_run = false;
  }

  // return stringstream to print current entries of
  // input_dict (can be send to logger)
  std::string Classy_input::print_entries_to_logger()
  {
    using namespace LogTags;
    std::ostringstream log_msg;
    log_msg << "Parameters passed to class: \n \n";
    for (auto item : input_dict)
    {
      pybind11::str key = pybind11::str(item.first);
      pybind11::str value = pybind11::str(item.second);

      log_msg << std::string(key) << " = " << std::string(value) << " \n";
    }
    return log_msg.str();
  }

  // clears all entries from input_dict
  void Classy_input::clear(){input_dict.attr("clear")();}

  // return input_dict
  pybind11::dict Classy_input::get_input_dict(){return input_dict;}

}
