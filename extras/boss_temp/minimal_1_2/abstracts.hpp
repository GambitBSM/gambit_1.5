#ifndef __abstracts_hpp__
#define __abstracts_hpp__

//Note that order is irrellevant
#include "backend_types/BOSSMinimalExample_1_2/abstract_Y.hpp"
#include "backend_types/BOSSMinimalExample_1_2/abstract_X.hpp"

#include "backend_types/BOSSMinimalExample_1_2/identification.hpp"
//Note that order is irrellevant
typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Abstract_X Abstract_X;
typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Abstract_Y Abstract_Y;
#include "backend_undefs.hpp"

#endif
