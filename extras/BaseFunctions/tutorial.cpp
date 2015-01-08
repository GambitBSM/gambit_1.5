#include "funktions.hpp"
// Compile with: 
//   g++ tutorial.cpp -lgsl -lblas -std=c++11

double plain_func(double x, std::string s, double& y)
{
    if ( s == "trig" )
    {
        return sin(x) * cos(y);
    }
    else
    {
        return x + y;
    }
}

class Object
{
public:
    Object() {};
    ~Object() {};

    double member_function(double x)
    {
        return x*x/2;
    }
};

DEF_FUNKTRAIT(T)  // Generate bucket for plain function without void pointer!

int main()
{
    {
        std::cout << "Creating functions from scratch" << std::endl;

        Funk::Funk x = Funk::var("x");
        x->help();
        std::cout << x->eval("x", 3.1415) << std::endl;

        Funk::Funk y = Funk::var("y");
        Funk::Funk f = x*y;
        f->help();
        std::cout << f->eval("x", 3.1415, "y", 3.1415) << std::endl;

        Funk::Funk g = Funk::one("x", "y")*2.71;
        g->help();
        std::cout << g->eval("x", 3.1415, "y", 3.1415) << std::endl;

        Funk::Funk h = Funk::zero("x", "y");
        g->help();
        std::cout << h->eval("x", 3.1415, "y", 3.1415) << std::endl;

        std::cout << std::endl;
    }

    {
        std::cout << "Creating functions from plain or member functions" << std::endl;

        Funk::Funk x = Funk::var("x");
        Funk::Funk y = Funk::var("y");
        Funk::Funk f = Funk::func(plain_func, x, "trig", y);
        f->help();
        std::cout << f->eval("x", 3.1415, "y", 3.1415) << std::endl;

        Object obj;
        Funk::Funk g = Funk::funcM(&obj, &Object::member_function, x);
        g->help();
        std::cout << g->eval("x", 3.1415) << std::endl;

        std::cout << std::endl;
    }

    {
        std::cout << "Mapping on plain functions" << std::endl;

        Funk::Funk x = Funk::var("x");
        Funk::Funk y = Funk::var("y");
        Funk::Funk f = x*y;

        double (*plain)(double&, double&) = f->plain<T>("x", "y");
        double a = 3.;
        double b = 4.;
        std::cout << plain(a, b) << std::endl;

        double (*plain_with_voidptr)(double, double, void*);
        void* ptr;
        std::tie(plain_with_voidptr, ptr) = f->plain("x", "y");
        std::cout << plain_with_voidptr(1, 1, ptr) << std::endl;

        std::cout << std::endl;
    }

    {
        std::cout << "Integration" << std::endl;

        Funk::Funk x = Funk::var("x");
        Funk::Funk y = Funk::var("y");
        Funk::Funk f = x*x*y*y;

        Funk::Funk g = f->gsl_integration("x", 0, y);
        std::cout << g->eval("y", 3) << std::endl;

        Funk::Funk h = f->gsl_integration("x", -1+y*y, 1-y*y);
        std::cout << h->gsl_integration("y", -1, 1)->eval() << std::endl;

        std::cout << std::endl;
    }

    {
        std::cout << "Mapping of variables" << std::endl;

        Funk::Funk p = Funk::var("p");
        Funk::Funk E = Funk::var("E");
        Funk::Funk v = Funk::var("v");
        Funk::Funk gamma = 1/sqrt(1-v*v);
        gamma->help();

        Funk::Funk Gamma = gamma->set("v", p/E);
        Gamma->help();

        double mass = 125.;
        Funk::Funk GAMMA = Gamma->set("p", sqrt(E*E-mass*mass));
        GAMMA->help();

        std::cout << std::endl;
    }
}
