#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

namespace nspace1
{

    namespace nspace2
    {

        } } 
        #include "backend_types/BOSSMinimalExample_1_0/abstract_X.hpp"
        #include "abstracts_typedefs.hpp"
        namespace nspace1 { namespace nspace2 { 
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
                Abstract_X* pointerCopy_GAMBIT();
                void pointerAssign_GAMBIT(Abstract_X* in);

            public:
                int& i_ref_GAMBIT();


                nspace1::nspace2::Abstract_X* operator_plus_GAMBIT(nspace1::nspace2::Abstract_X&);


            public:
                nspace1::nspace2::Abstract_X* return_ref_this_GAMBIT();

                nspace1::nspace2::Abstract_X* return_ptr_this_GAMBIT();

};

    }

}


namespace nspace3
{
    } 
    #include "backend_types/BOSSMinimalExample_1_0/abstract_Y.hpp"
    #include "abstracts_typedefs.hpp"
    namespace nspace3 { 
    class Y : public virtual Abstract_Y
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
    
        public:
            Abstract_Y* pointerCopy_GAMBIT();
            void pointerAssign_GAMBIT(Abstract_Y* in);

        public:
            nspace1::nspace2::Abstract_X& x_ref_GAMBIT();


        public:
            nspace1::nspace2::Abstract_X* get_x_GAMBIT();

            void set_x_GAMBIT(nspace1::nspace2::Abstract_X&);

            void set_x_ptr_GAMBIT(nspace1::nspace2::Abstract_X*);

};
}

#endif

