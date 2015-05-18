#include "../../modules/Elements/include/gambit/Elements/funktions.hpp"
#include <vector>
// Compile with: 
//   g++ tutorial.cpp -lgsl -lblas -std=c++11

// This does not actually crash, but was supposed to reproduce a segfault from
// the destructor of funktion objects.  Maybe it is useful later.

struct Dumb
{
    Funk::BoundFunk function;
};

Dumb makeDumb()
{
    Dumb dumb;
    auto x = Funk::var("x");
    auto y = Funk::var("y");
    auto f = x*y;
    dumb.function = f->bind("x", "y");
    return dumb;
}

int main()
{
    Dumb d = makeDumb();
}
