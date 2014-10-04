#ifndef __WRAPPERS_HPP__
#define __WRAPPERS_HPP__

#include "backend_types/BOSSMinimalExample_1_2/wrapper_X.hpp"
#include "backend_types/BOSSMinimalExample_1_2/wrapper_Y.hpp"
#include "backend_types/BOSSMinimalExample_1_2/identification.hpp"

typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::X wrapper_X;
typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::Y wrapper_Y;

#include "backend_undefs.hpp"

#endif /* __WRAPPERS_HPP__ */
