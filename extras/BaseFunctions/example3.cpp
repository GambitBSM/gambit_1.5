#include "funktions.hpp"

#include<iostream>
#include<memory>

template <typename R> class Funk;

template <typename R>
using ptr = std::shared_ptr<Funk<R>>;

template <typename R>
class Funk : public std::enable_shared_from_this<Funk<R>>
{
    public:
        virtual R value() = 0;

        ptr<R> blub()
        {
            return this->shared_from_this();
        }
};

template <typename R>
class FunkConst : public Funk<R>
{
    public:
        FunkConst(R x) : x(x) {};

        R value()
        {
            return x;
        }

    private:
        R x;
};
template <typename R>
ptr<R> con(R x)
{
    return std::shared_ptr<FunkConst<R>>(new FunkConst<R>(x));
}

template <typename R, typename X, typename Y>
class FunkSum : public Funk<R>
{
    public:
        FunkSum(ptr<X> & f1 , ptr<Y> & f2) : f1(f1), f2(f2) {};

        R value()
        {
            return f1->value() + f2->value();
        }

    private:
        ptr<X> f1;
        ptr<Y> f2;
};

template <typename X, typename Y>
ptr<X> sum(std::shared_ptr<Funk<X>> x, std::shared_ptr<Funk<Y>> y)
{
    return std::shared_ptr<FunkSum<X, X, Y>>(new FunkSum<X, X, Y>(x, y));
}

class Test
{
    public:
        double value(double r)
        {
            return r;
        }
};

template <typename T>
double process(T* obj, double (T::* f)(double))
{
    return (*obj.*f)(2);
}

int main()
{
    //ptr<double> c(new FunkConst<double>(3.));
    auto c = con(3.);
    auto s = sum(c, c);
    //ptr<double> s(new FunkSum<double,double,double>(c, c));
    //std::cout << s->value() << std::endl;
    Test test;
    std::cout << process(&test, &Test::value) << std::endl;
    //functor(test.value);
}
