#ifndef __GAMBIT_WRAPPER_X_HPP__
#define __GAMBIT_WRAPPER_X_HPP__

#include "../after_BOSS/abstract_X.hpp"
#include <iostream> // for testing

// Factory function pointers to be filled by dynamic loading
Abstract_X* (*Factory_X_1)() = NULL;
Abstract_X* (*Factory_X_2)(Abstract_T&) = NULL;


//
// Class X
//

class X_gambit
{
    private:
        bool member_variable;
    public:
        // Member variables:
        Abstract_X *BEptr;  
        T_gambit t;      

        // Member functions: 
        T_gambit getT()
        {
            return BEptr->getT();
        }

        void setT(T_gambit t_in)
        {
            BEptr->setT(*t_in.BEptr);
        }

        void refTest(T_gambit& t_in, int& i_in)
        {
            BEptr->refTest(*t_in.BEptr, i_in);
        }

        int**& testFunc(T_gambit* t1, T_gambit t2, int**& ipp, double d)
        {
            return BEptr->testFunc((*t1).BEptr, *t2.BEptr, ipp, d);
        }

        // Special member function to set member_variable: 
        void _set_member_variable(bool in) 
        { 
            member_variable = in; 
        }

        // FOR DEBUG: Special member funtion to print member_variable state
        void _print_member()
        {
            std::cout << " -- member_variable = " << member_variable << std::endl;
        }

        // Wrappers for original constructors: 
        X_gambit() :
            BEptr( Factory_X_1() ),
            t(&(BEptr->t_ref_GAMBIT())),
            member_variable(false)
        {
            std::cout << "(Default constructor of X_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
            t._set_member_variable(true);
        }

        X_gambit(T_gambit t_in) :
            BEptr( Factory_X_2(*t_in.BEptr) ),
            t(&(BEptr->t_ref_GAMBIT())),
            member_variable(false)
        {
            std::cout << "(Constructor (T_gambit) of X_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
            t._set_member_variable(true);
        }

        // Special pointer-based constructor: 
        X_gambit(Abstract_X* in) :
            BEptr(in),
            t(&(BEptr->t_ref_GAMBIT())),
            member_variable(false)
        {
            std::cout << "(Default (pointer-based) constructor of X_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
            t._set_member_variable(true);
        }

        // Copy constructor: 
        X_gambit(const X_gambit& in) :
            BEptr(in.BEptr->pointerCopy_GAMBIT()),
            t(&(BEptr->t_ref_GAMBIT())),
            member_variable(false)
        {
            std::cout << "(Copy constructor of X_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
            t._set_member_variable(true);
        }

        // Assignment operator: 
        X_gambit& operator=(const X_gambit& in)
        {
            if(member_variable==false)
            {
                std::cout << "(Destructor of X_gambit)" << std::endl;
                delete BEptr;
            }
        }

        // Destructor: 
        ~X_gambit()
        {
            if(member_variable==false)
            {
                std::cout << "(Destructor of T_gambit)" << std::endl;
                delete BEptr;
            }
        }
};

#endif /* __GAMBIT_WRAPPER_X_HPP__ */ 