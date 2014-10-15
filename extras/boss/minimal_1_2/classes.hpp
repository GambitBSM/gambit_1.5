#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

#include "backend_types/BOSSMinimalExample_1_2/abstract_X.hpp"
#include "abstracts_typedefs.hpp"
class X : public virtual Abstract_X 
{
    public:
        
        int i;

        X() : i(0) {}

        X(int i_in) : i(2*i_in) {}

        X& return_ref_this()
        {
            return *this;
        }

        X* return_ptr_this()
        {
            return this;
        }

        X operator+(X& x_rhs)
        {
            return X(i + x_rhs.i);
        }

    public:
        Abstract_X* pointerCopy_GAMBIT();
        void pointerAssign_GAMBIT(Abstract_X* in);

    public:
        int& i_ref_GAMBIT();


        Abstract_X* operator_plus_GAMBIT(Abstract_X&);


    public:
        Abstract_X* return_ref_this_GAMBIT();

        Abstract_X* return_ptr_this_GAMBIT();

};




#include "backend_types/BOSSMinimalExample_1_2/abstract_Y.hpp"
#include "abstracts_typedefs.hpp"
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

        void set_x_ptr(X* x_in)
        {
            x = *x_in;
        }

    public:
        Abstract_Y* pointerCopy_GAMBIT();
        void pointerAssign_GAMBIT(Abstract_Y* in);

    public:
        Abstract_X& x_ref_GAMBIT();


    public:
        Abstract_X* get_x_GAMBIT();

        void set_x_GAMBIT(Abstract_X&);

        void set_x_ptr_GAMBIT(Abstract_X*);

};


#endif

