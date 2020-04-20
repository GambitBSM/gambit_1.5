//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Declarations for the YAML options class.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Christoph Weniger
///          (c.weniger@uva.nl)
///  \date 2013 June 2013
///
///  \author Ben Farmer
///          (benjamin.farmer@monash.edu.au)
///  \date 2013 Dec
///
///  \author Gregory Martinez
///          (gregory.david.martinez@gmail.com)
///  \date 2014 Feb
///
///  \author Pat Scott
///          (patscott@physics.mcgill.ca)
///  \date 2014 Mar
///
///  *********************************************

#ifndef __yaml_options_hpp__
#define __yaml_options_hpp__

#include <vector>
#include <sstream>
#include <utility>

#include "gambit/Utils/util_types.hpp"
#include "gambit/Utils/standalone_error_handlers.hpp"
#include "gambit/Utils/yaml_variadic_functions.hpp"

namespace Gambit
{

  ///  A small wrapper object for 'options' nodes.
  ///  These can be extracted from the prior, observable/likelihood and rules sections of the
  ///  inifile, or set by hand in module standalone mode.
  class Options
  {

    public:

      /// Default constructor
      Options() {}

      /// Copy constructor
      Options(const YAML::Node &options) : options(options) {}

      /// Move constructor
      Options(YAML::Node &&options) : options(std::move(options)) {}

      /// Getters for key/value pairs (which is all the options node should contain)
      /// @{
      template <typename... args>
      bool hasKey(const args&... keys) const
      {
        return getVariadicNode(options, keys...);
      }

      template<typename TYPE, typename... args>
      TYPE getValue(const args&... keys) const
      {
        const YAML::Node node = getVariadicNode(options, keys...);
        TYPE result;
        if (not node)
        {
          std::ostringstream os;
          os << "No options entry for [" << stringifyVariadic(keys...) << "]\n Node contents:  " << options;
          utils_error().raise(LOCAL_INFO,os.str());
          result = TYPE();
        }
        else
        {
          try
          {
            result = getNode<TYPE>(node);
          }
          catch(YAML::Exception& e)
          {
            std::string nodestr;
            try
            {
              nodestr = node.as<std::string>();
            }
            catch(YAML::Exception& e)
            {
              nodestr = "<Couldn't even convert to string!>";
            }
            std::ostringstream os;
            os << "Error retrieving options entry for [" << stringifyVariadic(keys...)
               << "] as type " << typeid(TYPE).name() << " (template parameter: see below). String form of node value was: "
               << nodestr << std::endl
               << "YAML message follows: " << std::endl
               << e.what();
            utils_error().raise(LOCAL_INFO,os.str());
            result = TYPE();
          }
        }
        return result;
      }

      template<typename TYPE, typename... args>
      TYPE getValueOrDef(TYPE def, const args&... keys) const
      {
        const YAML::Node node = getVariadicNode(options, keys...);
        TYPE result;
        if (not node)
        {
          result = def;
        }
        else
        {
          result = getValue<TYPE>(keys...);
        }
        return result;
      }
      /// @}


      /// Basic setter, for adding extra options
      /// @{
      template<typename KEYTYPE, typename VALTYPE>
      void setValue(const KEYTYPE &key, const VALTYPE &val)
      {
         options[key] = val;
         return;
      }
      /// @}

      /// Retrieve values from key-value pairs in options node.
      /// Works for an arbitrary set of input keys (of any type), and returns
      /// all values as strings.
      template<typename... args>
      const std::vector<str> getNames(const args&... keys) const
      {
        std::vector<str> result;
        const YAML::Node node = getVariadicNode(options, keys...);

        if (node.IsMap())
        {
          for (auto it = node.begin(), end = node.end(); it != end; ++it)
          {
            result.push_back( it->first.as<str>() );
          }
        }

        return result;
      }

      /// Retrieve values from all key-value pairs in options node.
      /// Returns all keys as strings.
      const std::vector<str> getNames() const
      {
        std::vector<str> result;

        for (auto it = options.begin(), end = options.end(); it != end; ++it)
        {
          result.push_back( it->first.as<str>() );
        }

        return result;
      }

      /// Recursive options retrieval
      template<typename... args>
      const Options getOptions(const args&... keys) const
      {
        const YAML::Node node = getVariadicNode(options, keys...);
        if (node["options"])
        {
          return Options(node["options"]);
        }
        else
        {
          return Options(node);
        }
      }

      /// Retrieve raw YAML node
      template<typename... args>
      YAML::Node getNode(const args&... keys) const
      {
        const YAML::Node node = getVariadicNode(options, keys...);
        if (not node)
        {
          std::ostringstream os;
          os << "No options entry for [" << stringifyVariadic(keys...) << "]\n Node contents:  " << options;
          utils_error().raise(LOCAL_INFO,os.str());
        }
        return node;
      }

      /// Get YAML node from file
      template<typename... args>
      YAML::Node loadFromFile(const args&... keys) const
      {
        const YAML::Node node = getVariadicNode(options, keys...);
        if (not node)
        {
          std::ostringstream os;
          os << "No options entry for [" << stringifyVariadic(keys...) << "]\n Node contents:  " << options;
          utils_error().raise(LOCAL_INFO,os.str());
        }
        return YAML::LoadFile(node.as<std::string>().c_str());
      }

      /// Return begin and end of options
      YAML::const_iterator begin() const { return options.begin(); }
      YAML::const_iterator end() const { return options.end(); }

    private:

      YAML::Node options;

      /// Wrapper for reading the node for a given type. Default case does nothing.
      /// However in some instances we want to catch the yamlcpp exception and try
      /// to interpret it, e.g. scientific notation numbers as integers.
      template<class TYPE>
      TYPE getNode(const YAML::Node node) const { return node.as<TYPE>(); }

      template<class TYPE>
      TYPE safeIntegerTypeCast(const std::string& s) const {
        const double d = std::stod(s); // First savely convert to double and perform checks before casting to an integer type.
        if (d != std::floor(d))
        {
          std::ostringstream os;
          os << "Provided value " << d << " does not represent an integer.";
          utils_error().raise(LOCAL_INFO, os.str());
        }
        if (std::numeric_limits<TYPE>::max() < d or std::numeric_limits<TYPE>::min() > d)
        {
          std::ostringstream os;
          os << "Provided value " << d << " does not fit into the implemented integer type.";
          utils_error().raise(LOCAL_INFO, os.str());
        }
        return static_cast<TYPE>(d);
      }

  };

/// Allows to read scientific notation integer numbers. If the number does not
/// fit into the given type (here int) or is not an integer, this function will raise.
/// This exception is then caught by Options::getValue and handled.
template<>
inline int Options::getNode<int>(const YAML::Node node) const {
  try { return node.as<int>(); }
  catch (...) { return Options::safeIntegerTypeCast<int>(node.as<std::string>()); }
}

/// See int specialization.
template<>
inline unsigned int Options::getNode<unsigned int>(const YAML::Node node) const {
  try { return node.as<unsigned int>(); }
  catch (...) { return Options::safeIntegerTypeCast<unsigned int>(node.as<std::string>()); }
}

/// See int specialization.
template<>
inline long Options::getNode<long>(const YAML::Node node) const {
  try { return node.as<long>(); }
  catch (...) { return Options::safeIntegerTypeCast<long>(node.as<std::string>()); }
}

/// See int specialization.
template<>
inline unsigned long Options::getNode<unsigned long>(const YAML::Node node) const {
  try { return node.as<unsigned long>(); }
  catch (...) { return Options::safeIntegerTypeCast<unsigned long>(node.as<std::string>()); }
}

/// See int specialization.
template<>
inline long long Options::getNode<long long>(const YAML::Node node) const {
  try { return node.as<long long>(); }
  catch (...) { return Options::safeIntegerTypeCast<long long>(node.as<std::string>()); }
}

/// See int specialization.
template<>
inline unsigned long long Options::getNode<unsigned long long>(const YAML::Node node) const {
  try { return node.as<unsigned long long>(); }
  catch (...) { return Options::safeIntegerTypeCast<unsigned long long>(node.as<std::string>()); }
}

}

#endif //#ifndef __yaml_options_hpp__
