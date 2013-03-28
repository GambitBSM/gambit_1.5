/*
 * Proto backend rollcall
 * 
 * \author Pat Scott
 * \date 2013-03-27
 * 
 */

#ifndef __BACKEND_ROLLCALL_HPP__
#define __BACKEND_ROLLCALL_HPP__

#include "backend_libfirst.hpp"
#include "backend_libfortrancode.hpp"

#define BACKENDRENAME LibFirstCopy
  #include "backend_libfirst.hpp"
#undef BACKENDRENAME

#define BACKENDRENAME LibFortranCodeCopy
  #include "backend_libfortrancode.hpp"
#undef BACKENDRENAME

#endif /* __BACKEND_ROLLCALL_HPP__ */
