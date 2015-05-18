#include "../../modules/Elements/include/gambit/Elements/funktions.hpp"
#include <vector>
// Compile with: 
//   g++ tutorial.cpp -lgsl -lblas -std=c++11

// This does not actually crash, but was supposed to reproduce a segfault from
// the destructor of funktion objects.  Maybe it is useful later.

int main()
{
    auto x = Funk::var("x");
    auto y = Funk::var("y");
    auto f = x*y;
    auto F = f->bind("x", "y");
}
