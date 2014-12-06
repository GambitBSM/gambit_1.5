#include "funk.hpp"

double plain_f(double x, double y)
{
    return sin(x)*cos(x)*exp(y*y);
}

DEF_FUNKTRAIT(T)

int main()
{
    auto x = Funk::var("x");
    auto y = Funk::var("y");
    auto z = Funk::var("z");
    auto a = Funk::var("a");
    auto b = Funk::var("b");

    // Setting variables
    {
        auto f = exp(x) * log(y) * sin(z) + a + b;
        std::cout << f->eval("x", 1, "y", 2, "z", 3, "a", 4, "b", 5) << std::endl;
        std::cout << f->set("x", 1, "y", 2, "z", 3, "a", 4, "b", 5)->eval() << std::endl;
        std::cout << f->set("x", 1)->set("y", 2, "z", 3)->set("a", 4, "b", 5)->eval() << std::endl;
        std::cout << f->set("x", 1, "y", 2)->set("z", 3)->set("a", 4, "b", 5)->set("z", 5)->eval() << std::endl;
        std::cout << exp(1) * log(2) * sin(3) + 4 + 5 << std::endl;
    }

    // Substituting functions
    {
        auto f = exp(x) * log(y) * sin(z);
        auto g = f->help()->set("y", b)->help()->set("x", a)->help();
        std::cout << g->eval("z", 3, "a", 4, "b", 5) << std::endl;
        g = (exp(a) * log(b) * sin(z))->help();
        std::cout << g->eval("z", 3, "a", 4, "b", 5) << std::endl;
    }

    // Integration
    {
        double r = 2;
        double A = (x*x*y*y)->gsl_integration("y", -sqrt(r*r-x*x), sqrt(r*r-x*x))->gsl_integration("x", -r, r)->eval();
        std::cout << A << std::endl;
    }


    // Plain functions
    {
        auto f = Funk::func(plain_f, "x", "y");
        double r = 2;
        double A = f->gsl_integration("y", -sqrt(r*r-x*x), sqrt(r*r-x*x))->gsl_integration("x", -r, r)->eval();
        std::cout << A << std::endl;
        double (*plainF)(double);
        plainF = f->gsl_integration("y", -sqrt(a*a-x*x), sqrt(a*a-x*x))->gsl_integration("x", -a, a)->plain<T>("a");

        void *ptr;
        double(*plainF2)(double, void*);
        std::tie(plainF2, ptr) = f->gsl_integration("y", -sqrt(a*a-x*x), sqrt(a*a-x*x))->gsl_integration("x", -a, a)->plain("a");
        std::cout << plainF2(2, ptr) << std::endl;
        for (int i = 0; i < 1000; i++) {
            plainF2(1, ptr);
        }
    }
}
