#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

class X 
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
};




class Y
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
};
}

#endif

