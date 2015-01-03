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

template <typename... Args>
using PlainPtrs = std::pair<double(*)(Args...,void*), void*>;
template <typename... Args>
using PlainPtr = double(*)(Args...);

inline std::vector<double> linspace(double x0, double x1, unsigned int n)
{
    std::vector<double> ret;
    double dx = 0;
    if (n > 1)
        dx = (x1-x0)/(n-1);
        for (unsigned int i = 0; i<n; i++)
        {
            ret.push_back(x0 + i * dx);
        }
    return ret;
}

// Generate logarithmic 1-dim grid
inline std::vector<double> logspace(double x0, double x1, unsigned int n)
{
    std::vector<double> ret;
    double dx = 0;
    if (n > 1)
        dx = (x1-x0)/(n-1);
        for (unsigned int i = 0; i<n; i++)
        {
            ret.push_back(pow(10, x0 + i * dx));
        }
    return ret;
}


template <typename T>
std::vector<T> vec(std::vector<T> vector)
{
    return vector;
}
template <typename T, typename... Args>
std::vector<T> vec(std::vector<T> vector, T value, Args... args)
{
    vector.push_back(value);
    return vec(vector, args...);
}
// This function causes a (readable) compile-time error when T != U.
// In case types are convertable, they are converted.
template <typename T, typename U, typename... Args>
std::vector<T> vec(std::vector<T> vector, U value, Args... args)
{
    T value_converted = value;
    vector.push_back(value_converted);
    return vec(vector, args...);
}
template <typename T, typename... Args>
std::vector<T> vec(T value, Args... args)
{
    std::vector<T> vector;
    vector.push_back(value);
    vector = vec(vector, args...);
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

        PlainPtrs<double> plain(const char*);
        PlainPtrs<double,double> plain(const char*, const char*);
        PlainPtrs<double,double,double> plain(const char*, const char*, const char*);
        PlainPtrs<double,double,double,double> plain(const char*, const char*, const char*, const char*);

        template <typename T>
        PlainPtr<double> plain(const char*);
        template <typename T>
        PlainPtr<double,double> plain(const char*, const char*);
        template <typename T>
        PlainPtr<double,double,double> plain(const char*, const char*, const char*);
        template <typename T>
        PlainPtr<double,double,double,double> plain(const char*, const char*, const char*, const char*);

        virtual double value(const FunkArgs&) = 0;

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

class FunkPlain: public FunkBase
{
    public:
        FunkPlain(FunkPtr f, const char* arg1) : f(f), arg1(arg1) {};
        FunkPlain(FunkPtr f, const char* arg1, const char* arg2) : f(f), arg1(arg1), arg2(arg2) {};
        FunkPlain(FunkPtr f, const char* arg1, const char* arg2, const char* arg3) : f(f), arg1(arg1), arg2(arg2), arg3(arg3) {};
        FunkPlain(FunkPtr f, const char* arg1, const char* arg2, const char* arg3, const char* arg4) : f(f), arg1(arg1), arg2(arg2), arg3(arg3), arg4(arg4) {};

        static double plain2p(double x1, void* ptr)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1);
        }
        static double plain2p(double x1, double x2, void* ptr)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1, funkPtrPtr->arg2, x2);
        }
        static double plain3p(double x1, double x2, double x3, void* ptr)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1, funkPtrPtr->arg2, x2, funkPtrPtr->arg3, x3);
        }
        static double plain4p(double x1, double x2, double x3, double x4, void* ptr)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1, funkPtrPtr->arg2, x2, funkPtrPtr->arg3, x3, funkPtrPtr->arg4, x4);
        }

        template <typename T>
        static double plain1(double x1)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1);
        }
        template <typename T>
        static double plain2(double x1, double x2)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1, funkPtrPtr->arg2, x2);
        }
        template <typename T>
        static double plain3(double x1, double x2, double x3)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1, funkPtrPtr->arg2, x2, funkPtrPtr->arg3, x3);
        }
        template <typename T>
        static double plain4(double x1, double x2, double x3, double x4)
        {
            FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
            return funkPtrPtr->eval(funkPtrPtr->arg1, x1, funkPtrPtr->arg2, x2, funkPtrPtr->arg3, x3, funkPtrPtr->arg4, x4);
        }

        double value(const FunkArgs& args) 
        { 
            return f->value(args);
        }

    private:
        FunkPtr f;
        const char* arg1;
        const char* arg2;
        const char* arg3;
        const char* arg4;
};

PlainPtrs<double> FunkBase::plain(const char* arg1)
{
    void* ptr = new FunkPlain(shared_from_this(), arg1);
    return PlainPtrs<double>(&FunkPlain::plain2p, ptr);
}
PlainPtrs<double,double> FunkBase::plain(const char* arg1, const char* arg2)
{
    void* ptr = new FunkPlain(shared_from_this(), arg1, arg2);
    return PlainPtrs<double,double>(&FunkPlain::plain2p, ptr);
}
PlainPtrs<double,double,double> FunkBase::plain(const char* arg1, const char* arg2, const char* arg3)
{
    void* ptr = new FunkPlain(shared_from_this(), arg1, arg2, arg3);
    return PlainPtrs<double,double,double>(&FunkPlain::plain3p, ptr);
}
PlainPtrs<double,double,double,double> FunkBase::plain(const char* arg1, const char* arg2, const char* arg3, const char* arg4)
{
    void* ptr = new FunkPlain(shared_from_this(), arg1, arg2, arg3, arg4);
    return PlainPtrs<double,double,double,double>(&FunkPlain::plain4p, ptr);
}

template <typename T>
PlainPtr<double> FunkBase::plain(const char* arg1)
{
    T::ptr = new FunkPlain(shared_from_this(), arg1);
    return &FunkPlain::plain1<T>;
}
template <typename T>
PlainPtr<double,double> FunkBase::plain(const char* arg1, const char* arg2)
{
    T::ptr = new FunkPlain(shared_from_this(), arg1, arg2);
    return &FunkPlain::plain2<T>;
}
template <typename T>
PlainPtr<double,double,double> FunkBase::plain(const char* arg1, const char* arg2, const char* arg3)
{
    T::ptr = new FunkPlain(shared_from_this(), arg1, arg2, arg3);
    return &FunkPlain::plain3<T>;
}
template <typename T>
PlainPtr<double,double,double,double> FunkBase::plain(const char* arg1, const char* arg2, const char* arg3, const char* arg4)
{
    T::ptr = new FunkPlain(shared_from_this(), arg1, arg2, arg3, arg4);
    return &FunkPlain::plain4<T>;
}

class FunkDerived: public FunkBase
{
    public:
        FunkDerived(FunkPtr f, FunkArgs argmap) : f(f), argmap(argmap) {};
        FunkDerived(FunkPtr f, FunkArgs argmap, std::map<const char*, FunkPtr> funmap) : f(f), argmap(argmap), funmap(funmap) {};

        double value(const FunkArgs & args)
        {
            // args dominate over myargs
            FunkArgs myargs = args;  // Make copy
            for ( auto it = funmap.begin(); it != funmap.end(); it ++ )
            {
                myargs[it->first] = (it->second)->value(myargs);
            }
            myargs.insert(argmap.begin(), argmap.end());
            return f->value(myargs);
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
        double value(const FunkArgs& args)                                                         \
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
        double value(const FunkArgs& args)                                                         \
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

#define MATH_OPERATION(MATH)                                                 \
class FunkMath_##MATH: public FunkBase                                       \
{                                                                            \
    public:                                                                  \
        FunkMath_##MATH(FunkPtr f) : f(f) {};                                \
        double value(const FunkArgs& args) { return MATH(f->value(args)); }  \
    private:                                                                 \
        FunkPtr f;                                                           \
};                                                                           \
FunkPtr MATH(FunkPtr f) { return FunkPtr(new FunkMath_##MATH(f)); }
MATH_OPERATION(log)
MATH_OPERATION(sqrt)
MATH_OPERATION(log10)
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
class FunkFunk: public FunkBase
{
    public:
        template <typename... Args>
        FunkFunk(double (*f)(funcargs...), Args... args)
        {
            ptr = f;
            argnames = vec(args...);
        }

        double value(const FunkArgs& arg)
        {
            return value_from_map(arg);
        }

    private:
        std::vector<const char*> argnames;
        FunkArgs argmap;
        double (*ptr)(funcargs...);

        // Return value according to map
        template <typename... Args>
        double value_from_map(const FunkArgs & argmap, Args&... args)
        {
            auto it = argmap.find(argnames[sizeof...(args)]);
            if (it != argmap.end())
            {
                double x = it->second;
                return value_from_map(argmap, args..., x);
            }
            else
            {
                std::cout << "ERROR. Missing parameter: " << argnames[sizeof...(args)] << std::endl;
                exit(-1);
            }
        }
        double value_from_map(const FunkArgs & argmap, funcargs&... args)
        {
            return (*ptr)(args...);
        }
};

template <typename... funcargs, typename... Args>
FunkPtr func(double (*f)(funcargs...), Args... args) {
    return FunkPtr(new FunkFunk<funcargs...>(f, args...));
}

class FunkVar: public FunkBase
{
    public:
        FunkVar(const char* var) : var(var) {};

        double value(const FunkArgs& args) { 
            auto it = args.find(var);
            if (it != args.end())
            {
                return args.at(var); 
            }
            else
            {
                std::cout << "ERROR. Missing parameter: " << var << std::endl;
                exit(-1);
            }
        }

    private:
        const char* var;
};

FunkPtr v(const char* var) { return FunkPtr(new FunkVar(var)); }

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

        double value(const FunkArgs& args)
        {
            auto it = args.find(arg);
            if (it != args.end())
            {
                return (this->*ptr)(args.at(arg));
            }
            else
            {
                std::cout << "ERROR. Missing parameter: " << arg << std::endl;
                exit(-1);
            }
        }

    private:
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
        
        double(FunkInterp::*ptr)(double);
        const char * arg;
        std::vector<double> Xgrid;
        std::vector<double> Ygrid;
        std::string mode;
};

FunkPtr interp(const char * arg, std::vector<double> x, std::vector<double> y) { return FunkPtr(new FunkInterp(arg, x, y)); }


double f(double x1, double x2)
{
    return x1 + x2*2;
}

double g(double x1, double x2, double x3, double x4)
{
    return pow(x1, 2) - pow(x2, 2) - pow(x3, 2) - pow(x4, 2);
}

double h(double x)
{
    return sqrt(x);
}

#define DEFINE_FUNK_TRAIT(C) class C { public: static void* ptr; }; void* C::ptr = 0;

DEFINE_FUNK_TRAIT(FUNK_F1)

int main()
{
    FunkPtr b = func(g, "x1", "x2", "x3", "x4")->set("x1", 3);

    FunkPtr n = v("x1") + 3;
    FunkPtr m = v("mass") + 3;

    std::cout << func(h, "x")->set("x", interp("x", vec(1., 2., 3.), vec(3., 4., 5.)))->eval("x", 2.121) << std::endl;
    std::cout << log10(func(h, "x")->set("x", interp("x", linspace(1, 3, 10), linspace(20, 40, 10))))->eval("x", 2.121) << std::endl;
    /*
    std::cout << "1" << std::endl;
    auto p = (pow(sin(v("x")),2) + pow(cos(v("y")),2)) * v("z");
    for (int i = 0; i < 1000000; i++)
    {
        p->eval("x", 0.5, "y", 0.5, "z", 0.5);
    }
    std::cout << "1" << std::endl;
    for (int i = 0; i < 1000000; i++)
    {
        ((pow(sin(v("x")),2) + pow(cos(v("y")),2)) * v("z"))->eval("x", 0.5, "y", 0.5, "z", 0.5);
    }
    std::cout << "1" << std::endl;
    */

    double(*p1)(double,double,void*);
    void* p2;

    std::tie(p1, p2) = (v("x")*v("y"))->plain("x", "y");
    std::cout << p1(3, 2, p2) << std::endl;

    //double(*p3)(double,double);
    auto p3 = (log(v("z"))*v("x")*v("y"))->plain<FUNK_F1>("x", "y", "z");
    std::cout << p3(4, 1, -1) << std::endl;

    auto f = log(v("z"));

    f->integrate("s", 3, 5);
}
