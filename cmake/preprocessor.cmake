# GAMBIT: Global and Modular BSM Inference Tool  
#************************************************
# \file                                          
#                                                
#  Cmake preprocessor variable configuration. 
#    
#************************************************
#                                                
#  Authors (add name and date if you modify):                                    
#                                                
#  \author Antje Putze
#          (antje.putze@lapth.cnrs.fr)              
#  \date 2015 Apr, Sep
#
#  \author Pat Scott
#          (p.scott@imperial.ac.uk)              
#  \date 2015 Apr
#                                               
#  \author Tomas Gonzalo
#          (t.e.gonzalo@fys.uio.no)
#  \date 2016 Sep
#
#************************************************

include(CheckIncludeFiles)

message("${Yellow}-- Updating GAMBIT with config data${ColourReset}")

# Check for STL regex
if(COMPILER_SUPPORTS_CXX11 OR COMPILER_SUPPORTS_CXX0X)
  check_include_files(regex.h HAVE_REGEX_H)
endif()

# Check for link.h
check_include_files(link.h HAVE_LINK_H)

# Check for Graphviz
find_program(GRAPHVIZ_FOUND dot)
if (GRAPHVIZ_FOUND)
  message("${BoldYellow}   Found graphviz.${ColourReset} Model and module function hierarchy plots will be enabled.")
else()
  message("${BoldRed}   Did not find graphviz. Model and module function hierarchy plots will not be produced.${ColourReset}")
endif()

# Define HAVE_GRAPHVIZ compiler option
include(CMakeDependentOption)
cmake_dependent_option(HAVE_GRAPHVIZ "create Graphviz files" ON "GRAPHVIZ_FOUND" OFF)

# Check for Mathematica
include(cmake/FindMathematica.cmake)
if(Mathematica_FOUND)
  message("${BoldYellow}   Found Mathematica")
  if(Mathematica_WSTP_FOUND)
    message("${BoldYellow}   Found Wolfram Symbolic Transfer Protocol. Mathematica based backends enabled")
    set(HAVE_MATHEMATICA 1)
    set(MATHEMATICA_WSTP_H "${Mathematica_WSTP_INCLUDE_DIR}/wstp.h")
  else()
    message("${BoldRed}  WSTP not found. Please make sure it is installed before attempting to use Mathematica backends")
    set(HAVE_MATHEMATICA 0)
  endif()
else()
  message("${BoldRed}   Mathematica not found. Backends using Mathematica will be disabled")
  set(HAVE_MATHEMATICA 0)
endif()
#cmake_dependent_option(HAVE_MATHEMATICA "Mathematica" ON "Mathematica_WSTP_FOUND" OFF)

# Configure the file
set(outhdr "${PROJECT_SOURCE_DIR}/cmake/include/gambit/cmake/cmake_variables.hpp")
configure_file("${PROJECT_SOURCE_DIR}/cmake/cmake_variables.hpp.in" ${outhdr})
message("${BoldYellow}   Configured ${outhdr}${ColourReset}")

message("${Yellow}-- Updating GAMBIT with config data - done.${ColourReset}")

