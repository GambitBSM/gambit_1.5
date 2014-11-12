#ifndef __CLASSES_HPP__
#define __CLASSES_HPP__

namespace nspace1
{

    namespace nspace2
    {

        } } 
        #include "backend_types/BOSSMinimalExample_1_0/abstract_X.hpp"
        #include "abstracttypedefs.hpp"
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
                Abstract_X* pointerCopy__BOSS();
                void pointerAssign__BOSS(Abstract_X* in);

            public:
                int& i_ref__BOSS();


                nspace1::nspace2::Abstract_X* operator_plus__BOSS(nspace1::nspace2::Abstract_X&);


            public:
                nspace1::nspace2::Abstract_X* return_ref_this__BOSS();

                nspace1::nspace2::Abstract_X* return_ptr_this__BOSS();

};

    }

}


namespace nspace3
{
    } 
    #include "backend_types/BOSSMinimalExample_1_0/abstract_Y.hpp"
    #include "abstracttypedefs.hpp"
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
            Abstract_Y* pointerCopy__BOSS();
            void pointerAssign__BOSS(Abstract_Y* in);

        public:
            nspace1::nspace2::Abstract_X& x_ref__BOSS();


        public:
            nspace1::nspace2::Abstract_X* get_x__BOSS();

            void set_x__BOSS(nspace1::nspace2::Abstract_X&);

            void set_x_ptr__BOSS(nspace1::nspace2::Abstract_X*);

};
}

#endif

