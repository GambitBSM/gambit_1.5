#ifndef __WRAPPERS_HPP__
#define __WRAPPERS_HPP__

#include "backend_types/BOSSMinimalExample_1_0/wrapper_X.hpp"
#include "backend_types/BOSSMinimalExample_1_0/wrapper_Y.hpp"
#include "backend_types/BOSSMinimalExample_1_0/identification.hpp"

namespace nspace1{ namespace nspace2 { typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace1::nspace2::X wrapper_X; }}
namespace nspace3{ typedef CAT_3(BACKENDNAME,_,SAFE_VERSION)::nspace3::Y wrapper_Y; }

#include "backend_undefs.hpp"

#endif /* __WRAPPERS_HPP__ */
