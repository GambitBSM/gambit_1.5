#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__


namespace NamespaceForY
{
    class Y;
}


namespace NamespaceForX
{

    class X 
    {
        public:
            
            int i;
            NamespaceForY::Y* yptr;
            double* dptr;


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

            NamespaceForY::Y* get_yptr()
            {
                return yptr;
            }

            void set_yptr(NamespaceForY::Y* yptr_in)
            {
                yptr = yptr_in;
            }

            X operator+(X& x_rhs)
            {
                return X(i + x_rhs.i);
            }
    };
}



namespace NamespaceForY
{
    class Y
    {
        public:

            NamespaceForX::X x;

            Y() {}

            Y(NamespaceForX::X x_in) : x(x_in)
            {
                x = x_in;
            }

            NamespaceForX::X get_x()
            {
                return x;
            }

            void set_x(NamespaceForX::X& x_in)
            {
                x = x_in;
            }

            void set_x_ptr(NamespaceForX::X* x_in)
            {
                x = *x_in;
            }
    };
}

#endif

