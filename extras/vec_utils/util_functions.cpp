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
  
#include "util_functions.hpp"
  
namespace Gambit
{
  
  namespace Utils
  {
    
    /// Function to help static initialisation of our const data member vectors.
    /// Returns a copy of the vector with the string argument appended.
    std::vector<str> vecappend(const std::vector<str>& basevec, const str& newentry)
    {
      std::vector<str> newvec(basevec);
      newvec.push_back(newentry);
      return newvec;
    }
      
    /// Similar to vecappend(); joins two vectors and returns the result
    std::vector<str> vecjoin(const std::vector<str>& bv1, 
                             const std::vector<str>& bv2) 
    {
      std::vector<str> newvec;
      newvec.reserve( bv1.size() + bv2.size() );
      newvec.insert( newvec.end(), bv1.begin(), bv1.end() );
      newvec.insert( newvec.end(), bv2.begin(), bv2.end() );
      return newvec;
    }
        
    /// As per vecjoin() but joins three vectors and returns the result.
    std::vector<str> vecjoin3(const std::vector<str>& bv1, 
                              const std::vector<str>& bv2,
                              const std::vector<str>& bv3) 
    {
      std::vector<str> newvec;
      newvec.reserve( bv1.size() + bv2.size() + bv3.size() );
      newvec.insert( newvec.end(), bv1.begin(), bv1.end() );
      newvec.insert( newvec.end(), bv2.begin(), bv2.end() );
      newvec.insert( newvec.end(), bv3.begin(), bv3.end() );
      return newvec;
    }
  
  }

}
