//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  GAMBIT executable.
///
///  *********************************************
///
///  Authors:
///
///  \author Pat Scott
///          p.scott@imperial.ac.uk
///  \date 2019 June
///
///  *********************************************

#include <cstdlib>
#include <string>
#include <sstream>

#include "gambit/cmake/cmake_variables.hpp"

// Initializer; runs as soon as this library is loaded.
__attribute__((constructor))
static void initializer()
{
  std::ostringstream outstream;
  outstream << "\n\x1b[1;33mGAMBIT " << GAMBIT_VERSION_MAJOR << "." << GAMBIT_VERSION_MINOR << "." << GAMBIT_VERSION_REVISION;
  std::string patch(GAMBIT_VERSION_PATCH);
  if (patch != "") outstream << "-" << patch;
  outstream << "\nhttp://gambit.hepforge.org\n\x1b[0m";
  printf("%s", outstream.str().c_str());
  #ifndef EXCLUDE_RESTFRAMES
    const char* oldenv = getenv("CPLUS_INCLUDE_PATH");
    std::string existing(oldenv ? oldenv : "");
    std::string newenv = std::string(RESTFRAMES_INCLUDE) + ":" + existing;
    setenv("CPLUS_INCLUDE_PATH", newenv.c_str(), 1);
  #endif
}