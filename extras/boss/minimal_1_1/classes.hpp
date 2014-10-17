#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

#include "backend_types/BOSSMinimalExample_1_1/abstract_X.hpp"
#include "abstracttypedefs.hpp"
class X : public virtual Abstract_X 
{
    public:
        
        int i;

        X() : i(0) {}

        X(int i_in) : i(i_in) {}

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
        Abstract_X* pointerCopy__BOSS();
        void pointerAssign__BOSS(Abstract_X* in);

    public:
        int& i_ref__BOSS();


        Abstract_X* operator_plus__BOSS(Abstract_X&);


    public:
        Abstract_X* return_ref_this__BOSS();

        Abstract_X* return_ptr_this__BOSS();

};




#include "backend_types/BOSSMinimalExample_1_1/abstract_Y.hpp"
#include "abstracttypedefs.hpp"
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
        Abstract_Y* pointerCopy__BOSS();
        void pointerAssign__BOSS(Abstract_Y* in);

    public:
        Abstract_X& x_ref__BOSS();


    public:
        Abstract_X* get_x__BOSS();

        void set_x__BOSS(Abstract_X&);

        void set_x_ptr__BOSS(Abstract_X*);

};


#endif

