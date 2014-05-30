#ifndef __GAMBIT_WRAPPER_T_HPP__
#define __GAMBIT_WRAPPER_T_HPP__

#include "../after_BOSS/abstract_T.hpp"
#include <iostream> // for testing

// Factory function pointers to be filled by dynamic loading
Abstract_T* (*Factory_T_1)() = NULL;
Abstract_T* (*Factory_T_2)(int, double) = NULL;


//
// Class T
//

class T_gambit
// class T_gambit : public U_gambit
{
    private:
        bool member_variable;
    public:
        // Member variables:
        Abstract_T *BEptr;  
        int &i;
        double &d;
        
        // Member functions:
        void printMe()
        {
            BEptr->printMe();
        }
        
        // -- Inherited: from class U
        void memberFunc()
        {
            BEptr->memberFunc();
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
        T_gambit() : 
            BEptr(Factory_T_1()),
            i(BEptr->i_ref_GAMBIT()), 
            d(BEptr->d_ref_GAMBIT()), 
            member_variable(false)
        {
            std::cout << "(Default constructor of T_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
        }

        T_gambit(int i_in, double d_in) :
            BEptr( Factory_T_2(i_in, d_in) ),
            i(BEptr->i_ref_GAMBIT()),
            d(BEptr->d_ref_GAMBIT()),
            member_variable(false)
        {
            std::cout << "(Constructor (int,double) of T_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
        }

        // Special pointer-based constructor: 
        T_gambit(Abstract_T* in) :
            BEptr(in),
            i(BEptr->i_ref_GAMBIT()),
            d(BEptr->d_ref_GAMBIT()),
            member_variable(false)
        {
            std::cout << "(Default (pointer-based) constructor of T_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
        }

        // Copy constructor: 
        T_gambit(const T_gambit& in) :
            BEptr(in.BEptr->pointerCopy_GAMBIT()),
            i(BEptr->i_ref_GAMBIT()),
            d(BEptr->d_ref_GAMBIT()),
            member_variable(false)
        {
            std::cout << "(Copy constructor of T_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
        }

        // Assignment operator: 
        T_gambit& operator=(const T_gambit& in)
        {
            if (this != &in) { BEptr->pointerAssign_GAMBIT(in.BEptr); }
        }

        // Destructor:
        virtual ~T_gambit()
        {
            if(member_variable==false)
            {
                std::cout << "(Destructor of T_gambit)" << std::endl;
                delete BEptr;
            }
        }
};

#endif /* __GAMBIT_WRAPPER_T_HPP__ */ 