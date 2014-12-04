#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <functional>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>

#include "boost/shared_ptr.hpp"
#include "boost/enable_shared_from_this.hpp"

using boost::shared_ptr;
using boost::enable_shared_from_this;

class FunkBase;

typedef shared_ptr<FunkBase> FunkPtr;
typedef std::map<const char*, double> FunkArgs;



template <typename T>
std::vector<T> initVector(std::vector<T> vector)
{
    return vector;
}
template <typename T, typename... Args>
std::vector<T> initVector(std::vector<T> vector, T value, Args... args)
{
    vector.push_back(value);
    return initVector(vector, args...);
}
// This function causes a (readable) compile-time error when T != U.
// In case types are convertable, they are converted.
template <typename T, typename U, typename... Args>
std::vector<T> initVector(std::vector<T> vector, U value, Args... args)
{
    T value_converted = value;
    vector.push_back(value_converted);
    return initVector(vector, args...);
}
template <typename T, typename... Args>
std::vector<T> initVector(T value, Args... args)
{
    std::vector<T> vector;
    vector.push_back(value);
    vector = initVector(vector, args...);
    return vector;
}




class FunkBase: public enable_shared_from_this<FunkBase>
{
    public:
        FunkBase() {}
        ~FunkBase() {}

        //double operator() (const Farg &args) { return this->value(args); }

        template <typename... Args> FunkPtr set (Args... args);
        template <typename... Args> double eval (Args... args);

        virtual double value(FunkArgs&) = 0;

    private:
        // Read argument list.
        FunkArgs tmp_argmap;
        std::map<const char*, FunkPtr> tmp_funmap;
        template <typename... Args>
        void digest_arguments(const char* arg, double y, Args... args)
        {
            tmp_argmap[arg] = y;
            digest_arguments(args...);
        }
        template <typename... Args>
        void digest_arguments(const char* arg, FunkPtr y, Args... args)
        {
            tmp_funmap[arg] = y;
            digest_arguments(args...);
        }
        void digest_arguments() {};
};

class FunkDerived: public FunkBase
{
    public:
        FunkDerived(FunkPtr f, FunkArgs argmap) : f(f), argmap(argmap) {};
        FunkDerived(FunkPtr f, FunkArgs argmap, std::map<const char*, FunkPtr> funmap) : f(f), argmap(argmap), funmap(funmap) {};

        double value(FunkArgs & args)
        {
            // args dominate over myargs
            for ( auto it = funmap.begin(); it != funmap.end(); it ++ )
            {
                args[it->first] = (it->second)->value(args);
            }
            args.insert(argmap.begin(), argmap.end());
            return f->value(args);
        }

    private:
        FunkPtr f;
        FunkArgs argmap;
        std::map<const char*, FunkPtr> funmap;
};

#define MATH_OPERATION(OPERATION, SYMBOL)                                                          \
class FunkMath_##OPERATION: public FunkBase                                                        \
{                                                                                                  \
    public:                                                                                        \
        FunkMath_##OPERATION(FunkPtr f1, FunkPtr f2) : f1(f1), f2(f2), mode(0) {};                 \
        FunkMath_##OPERATION(FunkPtr f,  double x) : f1(f) , x(x), mode(1) {};                     \
        FunkMath_##OPERATION(double x,  FunkPtr f) : f1(f),  x(x), mode(2) {};                     \
        double value(FunkArgs& args)                                                               \
        {                                                                                          \
            if ( mode == 0 ) return f1->value(args) SYMBOL f2->value(args);                        \
            else if( mode == 1) return f1->value(args) SYMBOL x;                                   \
            else return x SYMBOL f1->value(args);                                                  \
        }                                                                                          \
    private:                                                                                       \
        FunkPtr f1, f2;                                                                            \
        double x;                                                                                  \
        int mode;                                                                                  \
};                                                                                                 \
FunkPtr operator SYMBOL (FunkPtr f1, FunkPtr f2) { return FunkPtr(new FunkMath_##OPERATION(f1, f2)); }  \
FunkPtr operator SYMBOL (double  x1, FunkPtr f2) { return FunkPtr(new FunkMath_##OPERATION(x1, f2)); }  \
FunkPtr operator SYMBOL (FunkPtr f1,  double x2) { return FunkPtr(new FunkMath_##OPERATION(f1, x2)); }
MATH_OPERATION(Sum,+)
MATH_OPERATION(Mul,*)
MATH_OPERATION(Div,/)
MATH_OPERATION(Dif,-)
#undef MATH_OPERATION

#define MATH_OPERATION(MATH)                                                                       \
class FunkMath_##MATH: public FunkBase                                                             \
{                                                                                                  \
    public:                                                                                        \
        FunkMath_##MATH(FunkPtr f1, FunkPtr f2) : f1(f1), f2(f2), mode(0) {};                      \
        FunkMath_##MATH(FunkPtr f,  double x) : f1(f) , x(x), mode(1) {};                          \
        FunkMath_##MATH(double x,  FunkPtr f) : f1(f),  x(x), mode(2) {};                          \
        double value(FunkArgs& args)                                                               \
        {                                                                                          \
            if ( mode == 0 ) return MATH(f1->value(args), f2->value(args));                        \
            else if( mode == 1) return MATH(f1->value(args), x);                                   \
            else return MATH(x, f1->value(args));                                                  \
        }                                                                                          \
    private:                                                                                       \
        FunkPtr f1, f2;                                                                            \
        double x;                                                                                  \
        int mode;                                                                                  \
};                                                                                                 \
FunkPtr MATH (FunkPtr f1, FunkPtr f2) { return FunkPtr(new FunkMath_##MATH(f1, f2)); }  \
FunkPtr MATH (double  x1, FunkPtr f2) { return FunkPtr(new FunkMath_##MATH(x1, f2)); }  \
FunkPtr MATH (FunkPtr f1,  double x2) { return FunkPtr(new FunkMath_##MATH(f1, x2)); }
MATH_OPERATION(pow)
#undef MATH_OPERATION

#define MATH_OPERATION(MATH)                                           \
class FunkMath_##MATH: public FunkBase                                 \
{                                                                      \
    public:                                                            \
        FunkMath_##MATH(FunkPtr f) : f(f) {};                          \
        double value(FunkArgs& args) { return MATH(f->value(args)); }  \
    private:                                                           \
        FunkPtr f;                                                     \
};                                                                     \
FunkPtr MATH(FunkPtr f) { return FunkPtr(new FunkMath_##MATH(f)); }
MATH_OPERATION(log)
MATH_OPERATION(sin)
MATH_OPERATION(cos)
MATH_OPERATION(tan)
MATH_OPERATION(exp)
#undef MATH_OPERATION


template <typename... Args> FunkPtr FunkBase::set (Args... args)
{
    tmp_argmap.clear();
    tmp_funmap.clear();
    digest_arguments(args...);
    return FunkPtr(new FunkDerived(shared_from_this(), tmp_argmap, tmp_funmap));
}
template <typename... Args> double FunkBase::eval (Args... args)
{
    tmp_argmap.clear();
    tmp_funmap.clear();
    digest_arguments(args...);
    return value(tmp_argmap);
}


template <typename... funcargs>
class FunkPlain: public FunkBase
{
    public:
        template <typename... Args>
        FunkPlain(double (*f)(funcargs...), Args... args)
        {
            ptr = f;
            argnames = initVector(args...);
        }

        double value(FunkArgs& arg)
        {
            return value_from_map(arg);
        }

    private:
        std::vector<const char*> argnames;
        FunkArgs argmap;
        double (*ptr)(funcargs...);

        // Return value according to map
        template <typename... Args>
        double value_from_map(FunkArgs & argmap, Args&... args)
        {
            auto it = argmap.find(argnames[sizeof...(args)]);
            if (it != argmap.end())
            {
                return value_from_map(argmap, args..., it->second);
            }
            else
            {
                std::cout << "ERROR. Missing parameter " << argnames[sizeof...(args)] << std::endl;
                exit(-1);
            }
        }
        double value_from_map(FunkArgs & argmap, funcargs&... args)
        {
            return (*ptr)(args...);
        }
};

template <typename... funcargs, typename... Args>
FunkPtr func(double (*f)(funcargs...), Args... args) {
    return FunkPtr(new FunkPlain<funcargs...>(f, args...));
}

class FunkVar: public FunkBase
{
    public:
        FunkVar(const char* var) : var(var) {};
        double value(FunkArgs& args) { return args[var]; }

    private:
        const char* var;
};

FunkPtr v(const char* var) { return FunkPtr(new FunkVar(var)); }




double f(double x1, double x2)
{
    return x1 + x2*2;
}

double g(double x1, double x2, double x3, double x4)
{
    return pow(x1, 2) - pow(x2, 2) - pow(x3, 2) - pow(x4, 2);
}

class FunkInterp : public FunkBase
{
    public:
        FunkInterp(const char * arg, std::vector<double> & Xgrid, std::vector<double> & Ygrid, std::string mode = "lin")
        {
            this->arg = arg;
            this->Xgrid = Xgrid;
            this->Ygrid = Ygrid;
            if ( mode == "lin" ) this->ptr = &FunkInterp::linearInterp;
            else if ( mode == "log" ) this->ptr = &FunkInterp::logInterp;
        };

        double value(FunkArgs& args)
        {
            return (this->*ptr)(args[arg]);
        }

        double logInterp(double x)
        {
            // Linear interpolation in log-log space
            if (x<Xgrid.front() or x>Xgrid.back()) return 0;
            int i = 0; for (; Xgrid[i] < x; i++) {};  // Find index
            double x0 = Xgrid[i-1];
            double x1 = Xgrid[i];
            double y0 = Ygrid[i-1];
            double y1 = Ygrid[i];
            return y0 * exp(log(y1/y0) * log(x/x0) / log(x1/x0));
        }

        double linearInterp(double x)
        {
            // Linear interpolation in lin-lin space
            if (x<Xgrid.front() or x>Xgrid.back()) return 0;
            int i = 0; for (; Xgrid[i] < x; i++) {};  // Find index
            double x0 = Xgrid[i-1];
            double x1 = Xgrid[i];
            double y0 = Ygrid[i-1];
            double y1 = Ygrid[i];
            return y0 + (x-x0)/(x1-x0)*(y1-y0);
        }
        
    private:
        double(FunkInterp::*ptr)(double);
        const char * arg;
        std::vector<double> Xgrid;
        std::vector<double> Ygrid;
        std::string mode;
};

FunkPtr interp(const char * arg, std::vector<double> x, std::vector<double> y) { return FunkPtr(new FunkInterp(arg, x, y)); }


int main()
{
    //auto abs = funkPlain(f, "lambda", "mass");
    //std::cout << abs("lambda", 4, "mass", 3) << std::endl;
    FunkPtr b = func(g, "x1", "x2", "x3", "x4")->set("x1", 3);

    FunkPtr n = v("x1") + 3;
    FunkPtr m = v("mass") + 3;

    cos(m)->eval("x1", 3, "x2", 4);

    std::cout << interp("x", initVector(1., 2., 3.), initVector(3., 4., 5.))->eval("x", 2.121) << std::endl;
}
