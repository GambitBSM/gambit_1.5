#include "funktions.hpp"

double plain_f(double x, double y)
{
    return sin(x)*cos(x)*exp(y*y);
}

double integrand(double s)
{
    const double phi = 1.;
    const double rho_s = 0.3;
    const double r_s = 20;
    const double Rsun = 8.5;
    double r_los = sqrt(pow(s, 2) + pow(Rsun, 2) - 2*s*Rsun*cos(phi));
    double rho_dm = rho_s*pow(r_s, 3)/r_los/pow(r_los+r_s, 2);
    return rho_dm*rho_dm;
}

double rho_func(double r, double r_s, double rho_s)
{
    return rho_s*pow(r_s, 3)/r/pow(r+r_s, 2);
}

double r_los_func(double s, double Rsun, double phi)
{
    return sqrt(pow(s, 2) + pow(Rsun, 2) - 2*s*Rsun*cos(phi));
}

double rho_func2(double r)
{
    const double rho_s = 0.3;
    const double r_s = 20;
    return rho_s*pow(r_s, 3)/r/pow(r+r_s, 2);
}

double r_los_func2(double s)
{
    const double Rsun = 8.5;
    const double phi = 1;
    return sqrt(pow(s, 2) + pow(Rsun, 2) - 2*s*Rsun*cos(phi));
}

double fort(int& i, double& x, double& y)
{
    return x*y + i;
}

DEF_FUNKTRAIT(T)

int main()
{
    auto x = Funk::var("x");
    auto y = Funk::var("y");
    auto z = Funk::var("z");
    auto a = Funk::var("a");
    auto b = Funk::var("b");

    auto Fort = Funk::func(fort, 3., y, z);
    std::cout << Fort->help()->eval("y", 3, "z", 4) << std::endl;

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
        std::cout << "Plain functions:" << std::endl;
        auto f = Funk::func(plain_f, x, y);
        double r = 2;
        double A = f->gsl_integration("y", -sqrt(r*r-x*x), sqrt(r*r-x*x))->gsl_integration("x", -r, r)->eval();
        std::cout << A << std::endl;
        double (*plainF)(double);
        plainF = f->gsl_integration("y", -sqrt(a*a-x*x), sqrt(a*a-x*x))->gsl_integration("x", -a, a)->plain<T>("a");
        (void)plainF;

        void *ptr;
        double(*plainF2)(double, void*);
        std::tie(plainF2, ptr) = f->gsl_integration("y", -sqrt(a*a-x*x), sqrt(a*a-x*x))->gsl_integration("x", -a, a)->plain("a");
        std::cout << plainF2(2, ptr) << std::endl;
        for (int i = 0; i < 1000; i++) {
            plainF2(1, ptr);
        }
    }

    // Line-of-sight integrals
    {
        {
            std::cout << "1: slow." << std::endl;

            double rho_s = 0.3;
            double r_s = 20;
            double Rsun = 8.5;

            auto r = Funk::var("r");
            auto s = Funk::var("s");
            auto phi = Funk::var("phi");

            auto rho_dm = rho_s*pow(r_s, 3)/r/pow(r+r_s, 2);
            auto r_los = sqrt(pow(s, 2) + pow(Rsun, 2) - 2*s*Rsun*cos(phi));

            auto integral = pow(rho_dm,2)->set("r", r_los)->set("phi", 1)->gsl_integration("s", 0, 1000);
            for (int i = 0; i < 10000; i++) {
                integral->eval();
            }
        }

        {
            std::cout << "2: medium." << std::endl;

            double rho_s = 0.3;
            double r_s = 20;
            double Rsun = 8.5;
            double phi = 1;

            auto rho_dm = Funk::func(rho_func, Funk::var("r"), Funk::var("r_s"), Funk::var("rho_s"));
            auto r_los = Funk::func(r_los_func, Funk::var("s"), Funk::var("Rsun"), Funk::var("phi"));
            auto integral = pow(rho_dm, 2)->set("r_s", r_s, "rho_s", rho_s, "r", r_los)->set("Rsun", Rsun, "phi", phi)->gsl_integration("s", 0, 1000);
            for (int i = 0; i < 10000; i++) {
                integral->eval();
            }
        }

        {
            std::cout << "3: medium." << std::endl;

            auto rho_dm = Funk::func(rho_func2, Funk::var("r"));
            auto r_los = Funk::func(r_los_func2, Funk::var("s"));
            auto integral = pow(rho_dm, 2)->set("r", r_los)->gsl_integration("s", 0, 1000);
            for (int i = 0; i < 10000; i++) {
                integral->eval();
            }
        }

        {
            std::cout << "4: fast." << std::endl;

            auto integral = Funk::func(integrand, Funk::var("s"))->gsl_integration("s", 0, 1000);
            for (int i = 0; i < 10000; i++) {
                integral->eval();
            }
        }
    }
}
