#include "../../modules/Elements/include/gambit/Elements/funktions.hpp"
// Compile with: 
//   g++ tutorial.cpp -lgsl -lblas -std=c++11


#define FUNCTION exp(x) * sin(x);

double function(double x)
{
    return FUNCTION;
}

int main()
{
    std::cout << "Start" << std::endl;
    {
        auto x = Funk::var("x");
        auto y = Funk::var("y");
        auto f = FUNCTION;
        double result;
        for ( int i = 0; i < 30000; i++)
        {
            auto b = f->gsl_integration("x", y, y+30)->bind("y");
            result = b->eval(0);
        }
        std::cout << result << std::endl;
    }

    {
        auto x = Funk::var("x");
        auto y = Funk::var("y");
        auto f = Funk::func(function, x);
        auto b = f->gsl_integration("x", y, y+30)->bind("y");
        double result;
        for ( int i = 0; i < 30000; i++)
        {
            result = b->eval(0);
        }
        std::cout << result << std::endl;
    }
}
