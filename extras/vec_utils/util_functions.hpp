//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions retired from GAMBIT on Nov 11 2014.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Pat Scott  
///          (patscott@physics.mcgill.ca)
///  \date 2013 Apr, July, Aug, Dec
///  \date 2014 Mar
///
///  \author Ben Farmer
///          (benjamin.farmer@monash.edu.au)
///  \date 2013 May, June, July
///
///  *********************************************


#ifndef __util_functions_hpp__
#define __util_functions_hpp__

#include <vector>
#include "util_types.hpp"

namespace Gambit
{

  namespace Utils
  {

    /// Function to help static initialisation of our const data member vectors.
    /// Returns a copy of the vector with the string argument appended.
    std::vector<str> vecappend(const std::vector<str>&, const str&);
   
    /// Similar to vecappend(); joins two vectors and returns the result
    std::vector<str> vecjoin(const std::vector<str>&, const std::vector<str>&);
      
    /// As per vecjoin() but joins three vectors and returns the result.
    std::vector<str> vecjoin3(const std::vector<str>&, 
                            const std::vector<str>&,
                            const std::vector<str>&);

  }

}
        
#endif //defined __util_functions_hpp__
