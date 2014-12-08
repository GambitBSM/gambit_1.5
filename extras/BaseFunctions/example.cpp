#include "funktions.hpp"
// Compile with: 
//   g++ example.cpp -lgsl -lblas

double f(double x, double &y, std::string s, const char * c, std::string& sr, int i)
{
    std::cout << s << " " << c << " " << sr << std::endl;
    return sin(x) * cos(y) + i;
}

DEF_FUNKTRAIT(T)  // Generate bucket for plain function without void pointer!

int main()
{
    auto x = Funk::var("x");
    auto y = Funk::var("y");

    // Wrap plain functions while passing arbitrary arguments for non-double
    // types
    auto F = Funk::func(f, x, y, "Pass", "it", "all!", 1);
    F->help();  // Print function arguments
    double result = F->eval("x", 1, "y", 1);
    std::cout << result << std::endl;

    // Getting back plain functions without void pointer...
    // (the pointer is carried by T, on which a static member function of
    // Funk::Funk is templated)
    double (*plain)(double, double) = F->plain<T>("x", "y");
    std::cout << plain(1, 1) << std::endl;

    // ...and with the usual void* pointer
    double (*plainV)(double, double, void*);
    void* ptr;
    std::tie(plainV, ptr) = F->plain("x", "y");
    std::cout << plainV(1, 1, ptr) << std::endl;

    // Do fancy integration
    auto G = sin(x) * cos(x);
    G->help();  // Print function arguments
    result = G->gsl_integration("x", 0, y)->eval("y", 3.1415/2);
    std::cout << "int sin(x)*cos(x) = " << result << std::endl;

    // 2dim integration with non-trivial boundaries (calculate pi!)
    double r = 1;
    double A = (x*x*y*y*0+1)->gsl_integration("y", -sqrt(r*r-x*x), sqrt(r*r-x*x))->gsl_integration("x", -r, r)->eval();
    std::cout << "2dim integral over unit sphere = " << A << std::endl;

    // Parameter mapping
    auto H = x*y->set("y", x*x);
    std::cout << "x**3 (with x=2) = " << H->eval("x", 2) << std::endl;
}
