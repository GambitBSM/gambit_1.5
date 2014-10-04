#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

#include "abstracts.hpp"
#include <iostream>

class X : public virtual Abstract_X 
{
    public:
        
        int i;

        X() : i(0) { std::cout << "Making an X with regular constructor." << std::endl; }

        X(int i_in) : i(i_in) { std::cout << "Making an X with int constructor." << std::endl; }

    public:
        Abstract_X* pointerCopy_GAMBIT() { return new X(*this); }
        void pointerAssign_GAMBIT(Abstract_X* in) { *this = *dynamic_cast<X*>(in); }

    public:
        int& i_ref_GAMBIT() { return i; }

};


class Y : public virtual Abstract_Y
{
    public:

        X x;

        Y() {}

        Y(X x_in) : x(x_in)
        {
            x = x_in;
        }

        X get_x()
        {
            return x;
        }

        void set_x(X& x_in)
        {
            x = x_in;
        }

    public:
        Abstract_Y* pointerCopy_GAMBIT() { return new Y(*this); }
        void pointerAssign_GAMBIT(Abstract_Y* in) { *this = *reinterpret_cast<Y*>(in); }

    public:
        Abstract_X& x_ref_GAMBIT() { return x; }


    public:
        Abstract_X* get_x_GAMBIT()
        {
            return new X(get_x());
        }

        void set_x_GAMBIT(Abstract_X& x_in)
        {
            set_x(reinterpret_cast< X& >(x_in));
        }

};


#endif
