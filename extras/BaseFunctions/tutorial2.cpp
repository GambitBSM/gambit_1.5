#include "funktions2.hpp"
// Compile with: 
//   g++ tutorial.cpp -lgsl -lblas -std=c++11

int main()
{
#define TEST1(EXPR, X, Y, Z) \
    { \
        Funk::Funk x = Funk::var("x");\
        Funk::Funk y = Funk::var("y");\
        Funk::Funk z = Funk::var("z");\
        std::cout << (EXPR)->bind("x", "z", "y")->eval(X, Z, Y) << " = ";\
    }\
    {\
        double x = X;\
        double y = Y;\
        double z = Z;\
        std::cout << EXPR << std::endl;\
    }
#define TEST2(EXPR, X, Y) \
    { \
        Funk::Funk x = Funk::var("x");\
        Funk::Funk y = Funk::var("y");\
        Funk::Funk z = Funk::var("z");\
        std::cout << (EXPR)->set("x", "z", "y", "x", "z", "y")->bind("x", "y")->eval(Y, X) << " = ";\
    }\
    {\
        double x = X;\
        double y = Y;\
        std::cout << EXPR << std::endl;\
    }
    {
        TEST1(x*y*z, 1, 2, 3)
            TEST1(sin(x)*y*z/(z*3+exp(4)), 1, 2, 3)
            TEST1(pow(x, y)*z, 2, 3, 4)
            TEST2(x+sin(x)/y, 2, 1.5)

            Funk::Funk x = Funk::var("x")->set_singularity("x", Funk::var("a"), Funk::var("b"))->help();
        std::cout << (x*x*pow(x, 3)/x/x/x)->help()->gsl_integration("x", 0, 3)->bind("a", "b")->eval(3, 4) << std::endl;
        auto xt = Funk::vec<double>(1,2,3,4);
        auto yt = Funk::vec<double>(2,3,4,5);
        std::cout << Funk::interp("x", xt, yt)->set("x", "y")->bind("y")->eval(3.5) << std::endl;

        double sigma = 1e-8;
        Funk::Funk a = Funk::var("a");
        Funk::Funk x1 = Funk::var("x");
        auto f = 1/sigma/sqrt(2*3.1415)*exp(-(x1-a)*(x1-a)/sigma/sigma/2)->set_singularity("x", a, Funk::cnst(sigma));
        double r = f->help()->gsl_integration("x", -1, 1)->bind("a")->eval(9e-1);
        std::cout << r << std::endl;
        auto g = (a*Funk::var("x"))->set_singularity("x", Funk::cnst(3), Funk::cnst(3))->set("a", "y");
    }

    {
        auto x = Funk::var("x");
        auto y = Funk::var("y");
        auto f = (x*sin(y))->set_singularity("x", 3, 0.3);
        f->set_singularity("x", 3, 0.3);
        f->set_singularity("x", 3, 0.3);
        f->set_singularity("x", 3, 0.3);
        f->gsl_integration("x", 3, 4)->set_singularity_factor(2.)->set_singularity_factor(2.);
    }
    {
        auto r = Funk::var("r");
        auto f = pow(r, -0.2);
        double result = f->gsl_integration("r", 0.0000, 10000)->bind()->eval();
        std::cout << result << std::endl;
        r->assert_args(Funk::vec<std::string>("r"));
    }
}
