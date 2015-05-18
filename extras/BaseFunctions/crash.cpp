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
    auto z = Funk::var("z");
    auto f = x*y;
    auto g = sin(x)*z;
    auto h = y+3*exp(z);
    auto F = f->bind("x", "y");
    auto G = g->bind("x", "z");
    auto H = h->bind("y", "z");

    std::vector<Funk::BoundFunk> v;
    std::vector<Funk::BoundFunk> w;
    for (int i=0; i<10000; i++)
    {
        auto I = f->bind("x", "y");
        auto J = g->bind("x", "z");
        v.push_back(I);
        w.push_back(J);
    }
}
