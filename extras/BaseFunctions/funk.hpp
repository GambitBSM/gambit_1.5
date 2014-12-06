/// ______           _    _   _                 
/// |  ___|         | |  | | (_)                
/// | |_ _   _ _ __ | | _| |_ _  ___  _ __  ___ 
/// |  _| | | | '_ \| |/ / __| |/ _ \| '_ \/ __|
/// | | | |_| | | | |   <| |_| | (_) | | | \__ \
/// \_|  \__,_|_| |_|_|\_\\__|_|\___/|_| |_|___/
///
///         Christoph Weniger, Dec 2014
///             <c.weniger@uva.nl>

#ifndef __FUNK_HPP__
#define __FUNK_HPP__
                                            
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
using boost::shared_ptr;
using boost::enable_shared_from_this;

//#include <memory>
//using std::shared_ptr;
//using std::enable_shared_from_this;

namespace Funk
{
    //
    // Type declarations
    //

    class FunkBase;
    typedef shared_ptr<FunkBase> FunkPtr;
    typedef std::vector<const char *> ArgsType;

    template <typename... Args>
    using PlainPtr = double(*)(Args...);

    template <typename... Args>
    using PlainPtrs = std::pair<double(*)(Args...,void*), void*>;


    //
    // Vector initialization from argument list
    // Usage: std::vector<T> v = vec<T>(v1, v2, v3, ...);
    //

    template <typename T>
    std::vector<T> vec(std::vector<T> vector)
    {
        return vector;
    }
    template <typename T, typename U, typename... Args>
    std::vector<T> vec(std::vector<T> vector, U value, Args... args)
    {
        T value_converted = value;  // Implicit conversion, if possible.
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


    //
    // Helper functions for internal calculations
    //

    ArgsType joinArgs(ArgsType args1, ArgsType args2)
    {
        args1.insert(args1.end(), args2.begin(), args2.end());
        std::set<const char*> argsset(args1.begin(), args1.end());
        args1.assign(argsset.begin(), argsset.end());
        return args1;
    }

    std::vector<int> getMap(ArgsType argIn, ArgsType argOut)
    {
        std::vector<int> map(argOut.size());
        for ( unsigned int i = 0; i < argOut.size(); i++ )
        {
            map[i] = std::find(argIn.begin(), argIn.end(), argOut[i]) - argIn.begin();
        }
        return map;
    }

    void applyMap(std::vector<double> & Xout, const std::vector<int> & map, const std::vector<double> & Xin, int n)
    {
        for ( int i = 0; i < n; i++ )
        {
            Xout[i] = Xin[map[i]];
        }
    }

    void applyInvMap(std::vector<double> & Xout, const std::vector<int> & map, const std::vector<double> & Xin, int n)
    {
        for ( int i = 0; i < n; i++ )
        {
            Xout[map[i]] = Xin[i];
        }
    }

    int eraseArg(ArgsType & args, const char* arg)
    {
        auto it = std::find(args.begin(), args.end(), arg);
        assert (it!=args.end());
        args.erase(it);
        return it - args.begin();
    }


    //
    // Central virtual base class
    //

    class FunkBase: public enable_shared_from_this<FunkBase>
    {
        public:
            FunkBase() {}
            ~FunkBase() {}

            // Standard handles
            template <typename... Args> FunkPtr set(Args... args);
            template <typename... Args> FunkPtr bind(Args... argss);
            template <typename... Args> double eval(Args... args);
            template <typename... Args> double get(Args... argss);

            // Convenience functions
            const std::vector<const char*> & getArgs() { return args; };
            void help();

            // Plain function generators (up to four arguments)
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

            // Return value
            virtual double value(const std::vector<double> &) = 0;

        protected:
            ArgsType args;  // Argument names

        private:
            // Internal structure required for bind
            std::vector<int> bind_map;
            std::vector<double> Xout;
            int nout;

            std::map<const char*, double> tmp_argmap;
            std::map<const char*, FunkPtr> tmp_funmap;
            std::vector<double> empty;

            // Add arguments from set(...) handler to tmp_argmap and
            // tmp_funmap.
            template <typename... Args> void digest_arguments(const char* arg, double y, Args... args);
            template <typename... Args> void digest_arguments(const char* arg, FunkPtr y, Args... args);
            void digest_arguments() {};
    };


    //
    // Derived class with (templated) static member functions as plain function
    // prototypes.
    //

    class FunkPlain: public FunkBase
    {
        public:
            FunkPlain(FunkPtr fin, const char* arg1) : f(fin->bind(arg1)) {}
            FunkPlain(FunkPtr fin, const char* arg1, const char* arg2) : f(fin->bind(arg1, arg2)) {}
            FunkPlain(FunkPtr fin, const char* arg1, const char* arg2, const char* arg3) : f(fin->bind(arg1, arg2, arg3)) {}
            FunkPlain(FunkPtr fin, const char* arg1, const char* arg2, const char* arg3, const char* arg4) : f(fin->bind(arg1, arg2, arg3, arg4)) {}

            static double plain1p(double x1, void* ptr)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
                return funkPtrPtr->f->get(x1);
            }
            static double plain2p(double x1, double x2, void* ptr)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
                return funkPtrPtr->f->get(x1, x2);
            }
            static double plain3p(double x1, double x2, double x3, void* ptr)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
                return funkPtrPtr->f->get(x1, x2, x3);
            }
            static double plain4p(double x1, double x2, double x3, double x4, void* ptr)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(ptr);
                return funkPtrPtr->f->get(x1, x2, x3, x4);
            }

            template <typename T>
            static double plain1(double x1)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
                return funkPtrPtr->f->get(x1);
            }
            template <typename T>
            static double plain2(double x1, double x2)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
                return funkPtrPtr->f->get(x1, x2);
            }
            template <typename T>
            static double plain3(double x1, double x2, double x3)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
                return funkPtrPtr->f->get(x1, x2, x3);
            }
            template <typename T>
            static double plain4(double x1, double x2, double x3, double x4)
            {
                FunkPlain * funkPtrPtr = static_cast<FunkPlain*>(T::ptr);
                return funkPtrPtr->f->get(x1, x2, x3, x4);
            }

            double value(const std::vector<double> & args) 
            { 
                assert ( 0 == 1 );  // This function should never be called
                return 0;
            }

        private:
            FunkPtr f;
            const char *arg1, *arg2, *arg3, *arg4;
    };


    //
    // Derived class that implements setting of parameters
    //

    class FunkDerived: public FunkBase
    {
        public:
            // Set parameter to fixed value
            FunkDerived(FunkPtr f, const char* arg, double x) : f(f), x(x), mode(0)
            {
                args = f->getArgs();
                XoutF.resize(args.size());
                i = eraseArg(args, arg);
                mapF = getMap(f->getArgs(), args);
                nF = XoutF.size();
            };
            // Map parameter to other function
            FunkDerived(FunkPtr f, const char* arg, FunkPtr g) : f(f), g(g), mode(1)
            {
                ArgsType argsF = f->getArgs();
                ArgsType argsG = g->getArgs();
                XoutF.resize(argsF.size());
                XoutG.resize(argsG.size());
                i = eraseArg(argsF, arg);
                args = joinArgs(argsG, args);
                mapF = getMap(argsF, args);
                mapG = getMap(args, argsG);
                nF = XoutF.size();
                nG = XoutG.size();
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
                    applyInvMap(XoutF, mapF, Xin, nF-1);
                    XoutF[i] = g->value(XoutG);
                    return f->value(XoutF);
                }
            }

        private:
            std::vector<double> XoutF, XoutG;
            std::vector<int> mapF, mapG;
            FunkPtr f, g;
            double x;
            int i, mode, nF, nG;
    };


    //
    // Derived class for the import of plain functions
    //

    template <typename... funcargs>
    class FunkFunc: public FunkBase
    {
        public:
            template <typename... Args>
            FunkFunc(double (*f)(funcargs...), Args... argss)
            {
                ptr = f;
                args = vec<const char*>(argss...);
            }

            double value(const std::vector<double> & X)
            {
                return args_from_vec(X);
            }

        private:
            double (*ptr)(funcargs...);

            // Return value according to map
            template <typename... Args>
            double args_from_vec(const std::vector<double> & X, Args... args)
            {
                return args_from_vec(X, args..., X[sizeof...(args)]);
            }
            double args_from_vec(const std::vector<double> & X, funcargs... args)
            {
                return (*ptr)(args...);
            }
    };

    template <typename... funcargs, typename... Args>
    FunkPtr func(double (*f)(funcargs...), Args... args) {
        return FunkPtr(new FunkFunc<funcargs...>(f, args...));
    }


    //
    // Derived class that implements simple linear variable
    //

    class FunkVar: public FunkBase
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
    FunkPtr var(const char* arg) { return FunkPtr(new FunkVar(arg)); }


    //
    // Definition of FunkBase member functions
    //

    template <typename... Args> FunkPtr FunkBase::set (Args... args)
    {
        tmp_argmap.clear();
        tmp_funmap.clear();
        digest_arguments(args...);
        FunkPtr f = shared_from_this();
        for ( auto it = tmp_argmap.begin(); it != tmp_argmap.end(); it++)
        {
            auto args = f->getArgs();
            if ( std::find(args.begin(), args.end(), it->first) != args.end() )
                f = FunkPtr(new FunkDerived(f, it->first, it->second));
            else
            {
                std::cout << "Funk: Ignoring \"" << it->first << "\" = " << it->second << std::endl;
            }
        }
        for ( auto it = tmp_funmap.begin(); it != tmp_funmap.end(); it++)
        {
            auto args = f->getArgs();
            if ( std::find(args.begin(), args.end(), it->first) != args.end() )
                f = FunkPtr(new FunkDerived(f, it->first, it->second));
            else
            {
                std::cout << "Funk: Ignoring \"" << it->first << "\" = function" << std::endl;
            }
        }
        return f;
    }

    template <typename... Args> FunkPtr FunkBase::bind(Args... argss)
    {
        ArgsType bind_args = vec<const char*>(argss...);
        assert (std::set<const char*>(args.begin(), args.end()) == std::set<const char*>(bind_args.begin(), bind_args.end()));
        bind_map = getMap(args, bind_args);
        Xout.resize(bind_args.size());
        nout = Xout.size();
        return shared_from_this();
    }

    template <typename... Args> double FunkBase::eval (Args... args)
    {
        if ( sizeof...(args) != 0 )
        {
            FunkPtr f = set(args...);
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

    template <typename... Args> double FunkBase::get(Args... argss)
    {
        assert (bind_map.size() == nout);
        applyMap(Xout, bind_map, vec<double>(argss...), nout);
        return this->value(Xout);
    }

    void FunkBase::help()
    {
        std::cout << "Arguments:";
        for ( auto it = args.begin(); it != args.end(); it++ )
        {
            std::cout << " \"" << *it << "\"";
        }
        std::cout << std::endl;
    }

    template <typename... Args> void FunkBase::digest_arguments(const char* arg, double y, Args... args)
    {
        tmp_argmap[arg] = y;
        digest_arguments(args...);
    }

    template <typename... Args> void FunkBase::digest_arguments(const char* arg, FunkPtr y, Args... args)
    {
        tmp_funmap[arg] = y;
        digest_arguments(args...);
    }

    PlainPtrs<double> FunkBase::plain(const char* arg1)
    {
        void* ptr = new FunkPlain(shared_from_this(), arg1);
        return PlainPtrs<double>(&FunkPlain::plain1p, ptr);
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


    //
    // Mathematical functions from cmath
    //

    // Unary minus sign
    class FunkMath_umin: public FunkBase                                                           
    {                                                                                                     
        public:                                                                                           
            FunkMath_umin(FunkPtr f) : f(f)                                                        
            {                                                                                             
                args = f->getArgs();                                                                      
            }                                                                                             
            double value(const std::vector<double> & X)                                                   
            {                                                                                             
                return -(f->value(X));                                                            
            }                                                                                            
        private:                                                                                         
            FunkPtr f;                                                                                    
    };                                                                                                    
    FunkPtr operator - (FunkPtr f) { return FunkPtr(new FunkMath_umin(f)); }

    // Unary operations
#define MATH_OPERATION(OPERATION)                                                                         \
    class FunkMath_##OPERATION: public FunkBase                                                           \
    {                                                                                                     \
        public:                                                                                           \
            FunkMath_##OPERATION(FunkPtr f) : f(f)                                                        \
            {                                                                                             \
                args = f->getArgs();                                                                      \
            }                                                                                             \
            double value(const std::vector<double> & X)                                                   \
            {                                                                                             \
                return OPERATION(f->value(X));                                                            \
            }                                                                                             \
        private:                                                                                          \
            FunkPtr f;                                                                                    \
    };                                                                                                    \
    FunkPtr OPERATION (FunkPtr f) { return FunkPtr(new FunkMath_##OPERATION(f)); }
    MATH_OPERATION(cos)
    MATH_OPERATION(sin)
    MATH_OPERATION(tan)
    MATH_OPERATION(acos)
    MATH_OPERATION(asin)
    MATH_OPERATION(atan)
    MATH_OPERATION(cosh)
    MATH_OPERATION(sinh)
    MATH_OPERATION(tanh)
    MATH_OPERATION(acosh)
    MATH_OPERATION(asinh)
    MATH_OPERATION(atanh)
    MATH_OPERATION(exp)
    MATH_OPERATION(log)
    MATH_OPERATION(log10)
    MATH_OPERATION(sqrt)
    MATH_OPERATION(fabs)
#undef MATH_OPERATION

    // Standard binary operations
#define MATH_OPERATION(OPERATION, SYMBOL)                                                                 \
    class FunkMath_##OPERATION: public FunkBase                                                           \
    {                                                                                                     \
        public:                                                                                           \
            FunkMath_##OPERATION(FunkPtr f1, FunkPtr f2) : f1(f1), f2(f2), mode(0)                        \
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
            FunkMath_##OPERATION(double x1, FunkPtr f2) : x1(x1), f2(f2), mode(1)                         \
            {                                                                                             \
                args = f2->getArgs();                                                                     \
            }                                                                                             \
            FunkMath_##OPERATION(FunkPtr f1, double x2) : x2(x2), f1(f1), mode(2)                         \
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
            FunkPtr f1, f2;                                                                               \
            int mode;                                                                                     \
    };                                                                                                    \
    FunkPtr operator SYMBOL (FunkPtr f1, FunkPtr f2) { return FunkPtr(new FunkMath_##OPERATION(f1, f2)); }\
    FunkPtr operator SYMBOL (double x, FunkPtr f) { return FunkPtr(new FunkMath_##OPERATION(x, f)); }     \
    FunkPtr operator SYMBOL (FunkPtr f, double x) { return FunkPtr(new FunkMath_##OPERATION(x, f)); }
    MATH_OPERATION(Sum,+)
    MATH_OPERATION(Mul,*)
    MATH_OPERATION(Div,/)
    MATH_OPERATION(Dif,-)
#undef MATH_OPERATION

    // More binary operations
#define MATH_OPERATION(OPERATION)                                                                         \
    class FunkMath_##OPERATION: public FunkBase                                                           \
    {                                                                                                     \
        public:                                                                                           \
            FunkMath_##OPERATION(FunkPtr f1, FunkPtr f2) : f1(f1), f2(f2), mode(0)                        \
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
            FunkMath_##OPERATION(double x1, FunkPtr f2) : x1(x1), f2(f2), mode(1)                         \
            {                                                                                             \
                args = f2->getArgs();                                                                     \
            }                                                                                             \
            FunkMath_##OPERATION(FunkPtr f1, double x2) : x2(x2), f1(f1), mode(2)                         \
            {                                                                                             \
                args = f1->getArgs();                                                                     \
            }                                                                                             \
            double value(const std::vector<double> & Xin)                                                 \
            {                                                                                             \
                if ( mode == 0 )                                                                          \
                {                                                                                         \
                    applyMap(Xout1, map1, Xin, n1);                                                       \
                    applyMap(Xout2, map2, Xin, n2);                                                       \
                    return OPERATION(f1->value(Xout1), f2->value(Xout2));                                 \
                }                                                                                         \
                if ( mode == 1 )                                                                          \
                {                                                                                         \
                    return OPERATION(x1, f2->value(Xin));                                                 \
                }                                                                                         \
                else                                                                                      \
                {                                                                                         \
                    return OPERATION(f1->value(Xin), x2);                                                 \
                }                                                                                         \
            }                                                                                             \
        private:                                                                                          \
            double x1, x2;                                                                                \
            int n1, n2;                                                                                   \
            std::vector<int> map1, map2;                                                                  \
            std::vector<double> Xout1, Xout2;                                                             \
            FunkPtr f1, f2;                                                                               \
            int mode;                                                                                     \
    };                                                                                                    \
    FunkPtr OPERATION (FunkPtr f1, FunkPtr f2) { return FunkPtr(new FunkMath_##OPERATION(f1, f2)); }      \
    FunkPtr OPERATION (double x, FunkPtr f) { return FunkPtr(new FunkMath_##OPERATION(x, f)); }           \
    FunkPtr OPERATION (FunkPtr f, double x) { return FunkPtr(new FunkMath_##OPERATION(f, x)); }
    MATH_OPERATION(pow)
    MATH_OPERATION(fmin)
    MATH_OPERATION(fmax)
#undef MATH_OPERATION

}

#define DEF_FUNKTRAIT(C) class C { public: static void* ptr; }; void* C::ptr = 0;

#endif  // __FUNK_HPP__
