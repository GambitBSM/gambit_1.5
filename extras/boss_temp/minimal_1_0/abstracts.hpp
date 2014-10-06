#ifndef __abstracts_hpp__
#define __abstracts_hpp__

//Note that order is irrellevant
#include "backend_types/BOSSMinimalExample_1_0/abstract_Y.hpp"
#include "backend_types/BOSSMinimalExample_1_0/abstract_X.hpp"

#include "backend_types/BOSSMinimalExample_1_0/identification.hpp"
//Note that order is irrellevant
namespace nspace1{ namespace nspace2{ typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace1::nspace2::Abstract_X Abstract_X; }}
namespace nspace3{ typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace3::Abstract_Y Abstract_Y; }
#include "backend_undefs.hpp"

#endif
