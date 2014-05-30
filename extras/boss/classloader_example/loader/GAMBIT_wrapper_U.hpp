#ifndef __GAMBIT_WRAPPER_U_HPP__
#define __GAMBIT_WRAPPER_U_HPP__

#include "../after_BOSS/abstract_U.hpp"
#include <iostream> // for testing

// Factory function pointers to be filled by dynamic loading
Abstract_U* (*Factory_U_1)() = NULL;


//
// Class U
//

class U_gambit
{
    private:
        bool member_variable;
    public:
        // Member variables:
        Abstract_U *BEptr;  
        
        // Member functions:
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
        U_gambit() : 
            BEptr(Factory_U_1()),
            member_variable(false)
        {
            std::cout << "(Default constructor of U_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
        }

        // Special pointer-based constructor: 
        U_gambit(Abstract_U* in) :
            BEptr(in),
            member_variable(false)
        {
            std::cout << "(Default (pointer-based) constructor of U_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
        }

        // Copy constructor: 
        U_gambit(const U_gambit& in) :
            BEptr(in.BEptr->pointerCopy_GAMBIT()),
            member_variable(false)
        {
            std::cout << "(Copy constructor of U_gambit)" << std::endl;
            // run _set_member_variable(true) on any member variables that are themselves wrapped classes
        }

        // Assignment operator: 
        U_gambit& operator=(const U_gambit& in)
        {
            if (this != &in) { BEptr->pointerAssign_GAMBIT(in.BEptr); }
        }

        // Destructor:
        virtual ~U_gambit()
        {
            if(member_variable==false)
            {
                std::cout << "(Destructor of U_gambit)" << std::endl;
                delete BEptr;
            }
        }
};

#endif /* __GAMBIT_WRAPPER_U_HPP__ */ 