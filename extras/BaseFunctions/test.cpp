#include <vector>
#include <string>
#include "funktions.hpp"
// Compile with: 
//   g++ tutorial.cpp -lgsl -lblas -std=c++11

struct test
{
    test(std::vector<std::string> bullshit, Funk::Funk w) : w(w)
    {
        w->help();
        w->set("v", 1);
    }
    Funk::Funk w;
};

int main()
{
    test t(Funk::vec<std::string>("1", "2", "3"), Funk::var("v"));
}
