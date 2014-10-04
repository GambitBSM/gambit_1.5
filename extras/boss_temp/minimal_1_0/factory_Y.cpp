#include "classes.hpp"
#include "wrappers.hpp"

namespace nspace3
{

  Abstract_Y* Factory_Y()
  {
    return new Y();
  }

  Abstract_Y* Factory_Y(nspace1::nspace2::wrapper_X& x_in)
  {
    return new Y(reinterpret_cast< nspace1::nspace2::X& >(*x_in.BEptr));
  }

}


#include "backend_undefs.hpp"
