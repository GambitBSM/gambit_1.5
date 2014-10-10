#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

namespace nspace1
{

    namespace nspace2
    {

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

    }

}


namespace nspace3
{
    class Y
    {
        public:

            nspace1::nspace2::X x;

            Y() {}

            Y(nspace1::nspace2::X x_in) : x(x_in)
            {
                x = x_in;
            }

            nspace1::nspace2::X get_x()
            {
                return x;
            }

            void set_x(nspace1::nspace2::X& x_in)
            {
                x = x_in;
            }

            void set_x_ptr(nspace1::nspace2::X* x_in)
            {
                x = *x_in;
            }
    };
}

#endif

