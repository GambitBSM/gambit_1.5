#include "classes.hpp"

namespace nspace1
{

  namespace nspace2
  {

    Abstract_X* Factory_X()
    {
      return new X();
    }

    Abstract_X* Factory_X(int i_in)
    {
      return new X(i_in);
    }

  }

}

