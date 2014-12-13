#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>

//#define NDEBUG
#include <assert.h>

#include "boost/shared_ptr.hpp"
#include "boost/enable_shared_from_this.hpp"

//#include <memory>
//using std::shared_ptr;
//using std::enable_shared_from_this;

#define DEF_FUNKTRAIT(C) class C { public: static void* ptr; }; void* C::ptr = 0;

// Extensions
#include <gsl/gsl_integration.h>

namespace Funk
{
    //
    // Type declarations
    //
    using boost::shared_ptr;
    using boost::enable_shared_from_this;

    class FunkBase;
    typedef shared_ptr<FunkBase> Funk;
    typedef std::vector<const char *> ArgsType;

    template <typename R>
    class FunkBase2;
    template <typename R>
    using Funk2 = shared_ptr<FunkBase2<R>>;


    //
    // Vector initialization from argument list
    // Usage: std::vector<T> v = vec<T>(v1, v2, v3, ...);
    //

    template <typename T>
    inline std::vector<T> vec(std::vector<T> vector)
    {
        return vector;
    }
    template <typename T, typename U, typename... Args>
    inline std::vector<T> vec(std::vector<T> vector, U value, Args... args)
    {
        T value_converted = value;  // Implicit conversion, if possible.
        vector.push_back(value_converted);
        return vec(vector, args...);
    }
    template <typename T, typename... Args>
    inline std::vector<T> vec(T value, Args... args)
    {
        std::vector<T> vector;
        vector.push_back(value);
        vector = vec(vector, args...);
        return vector;
    }


    //
    // Helper functions for internal calculations
    //

    inline ArgsType joinArgs(ArgsType args1, ArgsType args2)
    {
        args1.insert(args1.end(), args2.begin(), args2.end());
        std::set<const char*> argsset(args1.begin(), args1.end());
        args1.assign(argsset.begin(), argsset.end());
        return args1;
    }

    inline std::vector<int> getMap(ArgsType argIn, ArgsType argOut)
    {
        std::vector<int> map(argOut.size());
        for ( unsigned int i = 0; i < argOut.size(); i++ )
        {
            map[i] = std::find(argIn.begin(), argIn.end(), argOut[i]) - argIn.begin();
        }
        return map;
    }

    inline void applyMap(std::vector<double> & Xout, const std::vector<int> & map, const std::vector<double> & Xin, int n)
    {
        for ( int i = 0; i < n; i++ )
        {
            Xout[i] = Xin[map[i]];
        }
    }

    inline void applyInvMap(std::vector<double> & Xout, const std::vector<int> & map, const std::vector<double> & Xin, int n)
    {
        for ( int i = 0; i < n; i++ )
        {
            Xout[map[i]] = Xin[i];
        }
    }

    inline int eraseArg(ArgsType & args, const char* arg)
    {
        auto it = std::find(args.begin(), args.end(), arg);
        assert (it!=args.end());
        args.erase(it);
        return it - args.begin();
    }

    inline bool hasArg(ArgsType & args, const char* arg)
    {
        auto it = std::find(args.begin(), args.end(), arg);
        return it!=args.end();
    }


    //
    // Index lists (taken from stackoverflow)
    //

    // The structure that encapsulates index lists
    template <size_t... Is>
    struct index_list
    {
    };

    // Collects internal details for generating index ranges [MIN, MAX)
    namespace detail
    {
        // Declare primary template for index range builder
        template <size_t MIN, size_t N, size_t... Is>
        struct range_builder;

        // Base step
        template <size_t MIN, size_t... Is>
        struct range_builder<MIN, MIN, Is...>
        {
            typedef index_list<Is...> type;
        };

        // Induction step
        template <size_t MIN, size_t N, size_t... Is>
        struct range_builder : public range_builder<MIN, N - 1, N - 1, Is...>
        {
        };
    }

    // Meta-function that returns a [MIN, MAX) index range
    template<size_t MIN, size_t MAX>
    using index_range = typename detail::range_builder<MIN, MAX>::type;


    //
    // Central virtual base class
    //

    template <typename R>
    class FunkBase2: public enable_shared_from_this<FunkBase2<R>>
    {
        public:
            FunkBase2<R>() {}
            ~FunkBase2<R>() {}

            // Standard handles
            template <typename... Args> Funk set(Args... args);
            template <typename... Args> Funk bind(Args... argss);
            template <typename... Args> double eval(Args... args);
            template <typename... Args> double get(Args... argss);
            template <typename... Args> double operator() (Args... argss) { return this->eval(argss...); }
            std::vector<double> vector(const char*, const std::vector<double>&);

            // Extension handles
            // TODO: Implement
            // - tabularize
            template <typename... Args> Funk gsl_integration(Args... args);

            // Convenience functions
            const std::vector<const char*> & getArgs() { return args; };
            Funk help();

            // Return value
            virtual double value(const std::vector<double> &) = 0;

        protected:
            ArgsType args;  // Argument names

        private:
            // Internal structure required for bind
            std::vector<int> bind_map;
            std::vector<double> Xout;
            unsigned int nout;

            std::map<const char*, double> tmp_argmap;
            std::map<const char*, Funk> tmp_funmap;
            std::vector<double> empty;

            // Add arguments from set(...) handler to tmp_argmap and
            // tmp_funmap.
            template <typename... Args> void digest_arguments(const char* arg, double y, Args... args);
            template <typename... Args> void digest_arguments(const char* arg, Funk y, Args... args);
            void digest_arguments() {};
    };


    //
    // Derived class that implements setting of parameters
    //

    template <typename R>
    class FunkDerived: public FunkBase2<FunkDerived<R>>
    {
        public:
            // Set parameter to fixed value
            FunkDerived<R>(Funk2<R> f, const char* arg, double x) : f(f), x(x), mode(0)
            {
                //args = f->getArgs();
                //XoutF.resize(args.size());
                //i = eraseArg(args, arg);
                //mapF = getMap(f->getArgs(), args);
                //nF = XoutF.size();
            };
            // Map parameter to other function
            FunkDerived(Funk2<R> f, const char* arg, Funk2<R> g) : f(f), g(g), mode(1)
            {
                ArgsType argsF = f->getArgs();
                ArgsType argsG = g->getArgs();
                XoutF.resize(argsF.size());
                XoutG.resize(argsG.size());
                nF = XoutF.size();
                nG = XoutG.size();
                i = eraseArg(argsF, arg);
                //args = joinArgs(argsG, argsF);
                //mapF = getMap(args, f->getArgs());
                //mapG = getMap(args, argsG);
            };

            double value(const std::vector<double> & Xin)
            {
                if ( mode == 0 )
                {
                    applyInvMap(XoutF, mapF, Xin, nF-1);
                    XoutF[i] = x;
                    return f->value(XoutF);
                }
                else
                {
                    applyMap(XoutG, mapG, Xin, nG);
                    applyMap(XoutF, mapF, Xin, nF);
                    XoutF[i] = g->value(XoutG);
                    return f->value(XoutF);
                }
            }

        private:
            std::vector<double> XoutF, XoutG;
            std::vector<int> mapF, mapG;
            Funk f, g;
            double x;
            int i, mode, nF, nG;
    };


    //
    // Derived class that implements simple linear variable
    //

    class FunkVar: public FunkBase2<FunkVar>
    {
        public:
            FunkVar(const char* arg)
            {
                args = vec(arg);
            }

            double value(const std::vector<double> & X)
            {
                return X[0];
            }
    };
    inline Funk var(const char* arg) { return Funk(new FunkVar(arg)); }


    //
    // Definition of FunkBase member functions
    //

    template <typename... Args> inline Funk FunkBase2::set (Args... args)
    {
        tmp_argmap.clear();
        tmp_funmap.clear();
        digest_arguments(args...);
        Funk f = shared_from_this();
        for ( auto it = tmp_argmap.begin(); it != tmp_argmap.end(); it++)
        {
            auto args = f->getArgs();
            if ( std::find(args.begin(), args.end(), it->first) != args.end() )
                f = Funk(new FunkDerived(f, it->first, it->second));
            else
            {
                std::cout << "Funk: Ignoring \"" << it->first << "\" = " << it->second << std::endl;
            }
        }
        for ( auto it = tmp_funmap.begin(); it != tmp_funmap.end(); it++)
        {
            auto args = f->getArgs();
            if ( std::find(args.begin(), args.end(), it->first) != args.end() )
            {
                f = Funk(new FunkDerived(f, it->first, it->second));
            }
            else
            {
                std::cout << "Funk: Ignoring \"" << it->first << "\" = function" << std::endl;
            }
        }
        return f;
    }

    template <typename... Args> inline double FunkBase2::eval (Args... args)
    {
        if ( sizeof...(args) != 0 )
        {
            Funk f = set(args...);
            if ( f->args.size() == 0 )
            {
                return f->value(this->empty);
            }
            {
                std::cout << "Missing parameters." << std::endl;
                exit(-1);
            }
        }
        else
        {
            if ( this->args.size() == 0 )
            {
                return this->value(this->empty);
            }
            else
            {
                std::cout << "Missing parameters." << std::endl;
                exit(-1);
            }
        }
    }

    template <typename... Args> inline void FunkBase2::digest_arguments(const char* arg, double y, Args... args)
    {
        tmp_argmap[arg] = y;
        digest_arguments(args...);
    }

    template <typename... Args> inline void FunkBase2::digest_arguments(const char* arg, Funk y, Args... args)
    {
        tmp_funmap[arg] = y;
        digest_arguments(args...);
    }

    /*

    // Standard binary operations
#define MATH_OPERATION(OPERATION, SYMBOL)                                                                 \
    class FunkMath_##OPERATION: public FunkBase                                                           \
    {                                                                                                     \
        public:                                                                                           \
            FunkMath_##OPERATION(Funk f1, Funk f2) : f1(f1), f2(f2), mode(0)                        \
            {                                                                                             \
                ArgsType args1 = f1->getArgs();                                                           \
                ArgsType args2 = f2->getArgs();                                                           \
                args = joinArgs(args1, args2);                                                            \
                map1 = getMap(args, args1);                                                               \
                map2 = getMap(args, args2);                                                               \
                Xout1.resize(args1.size());                                                               \
                Xout2.resize(args2.size());                                                               \
                n1 = Xout1.size();                                                                        \
                n2 = Xout2.size();                                                                        \
            }                                                                                             \
            FunkMath_##OPERATION(double x1, Funk f2) : x1(x1), f2(f2), mode(1)                         \
            {                                                                                             \
                args = f2->getArgs();                                                                     \
            }                                                                                             \
            FunkMath_##OPERATION(Funk f1, double x2) : x2(x2), f1(f1), mode(2)                         \
            {                                                                                             \
                args = f1->getArgs();                                                                     \
            }                                                                                             \
            double value(const std::vector<double> & Xin)                                                 \
            {                                                                                             \
                if ( mode == 0 )                                                                          \
                {                                                                                         \
                    applyMap(Xout1, map1, Xin, n1);                                                       \
                    applyMap(Xout2, map2, Xin, n2);                                                       \
                    return f1->value(Xout1) SYMBOL f2->value(Xout2);                                      \
                }                                                                                         \
                if ( mode == 1 )                                                                          \
                {                                                                                         \
                    return x1 SYMBOL f2->value(Xin);                                                      \
                }                                                                                         \
                else                                                                                      \
                {                                                                                         \
                    return f1->value(Xin) SYMBOL x2;                                                      \
                }                                                                                         \
            }                                                                                             \
        private:                                                                                          \
            double x1, x2;                                                                                \
            int n1, n2;                                                                                   \
            std::vector<int> map1, map2;                                                                  \
            std::vector<double> Xout1, Xout2;                                                             \
            Funk f1, f2;                                                                               \
            int mode;                                                                                     \
    };                                                                                                    \
    inline Funk operator SYMBOL (Funk f1, Funk f2) { return Funk(new FunkMath_##OPERATION(f1, f2)); }\
    inline Funk operator SYMBOL (double x, Funk f) { return Funk(new FunkMath_##OPERATION(x, f)); }     \
    inline Funk operator SYMBOL (Funk f, double x) { return Funk(new FunkMath_##OPERATION(f, x)); }
    MATH_OPERATION(Sum,+)
    MATH_OPERATION(Mul,*)
    MATH_OPERATION(Div,/)
    MATH_OPERATION(Dif,-)
#undef MATH_OPERATION
*/
}
