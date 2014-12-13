#include "funktions.hpp"

struct Test
{
    double something(double r)
    {
        return r;
    }
};

int main()
{
    Test test;
    auto t = Funk::funcM(&test, &Test::something, Funk::var("r"));
    t->eval("r", 3);
}
