#ifndef __GAMBIT_WRAPPER_CLASSES_HPP__
#define __GAMBIT_WRAPPER_CLASSES_HPP__


#include "../after_BOSS/abstract_classes.hpp"
#include <iostream> // for testing

// Factory function pointers to be filled by dynamic loading
Abstract_T* (*Factory_T)() = NULL;
Abstract_X* (*Factory_X)() = NULL;
Abstract_Container<int>* (*Factory_Container_int)() = NULL;
Abstract_Container<Abstract_X>* (*Factory_Container_X)() = NULL;
Abstract_Container<Abstract_T>* (*Factory_Container_T)() = NULL;



//
// Class T
//

class T_gambit
{
    private:
        bool member_variable;
    public:
        Abstract_T *BEptr;  

        int &i;
        double &d;
        void printMe(){BEptr->printMe();};  

        T_gambit();
        T_gambit(Abstract_T*);        
        T_gambit(const T_gambit&);
        T_gambit& operator=(const T_gambit&);           
        virtual ~T_gambit();

        void _set_member(bool in) { member_variable = in; }
        void _print_member() { std::cout << " -- member_variable = " << member_variable << std::endl; }
};


// T_gambit class functions

T_gambit::T_gambit() : 
    BEptr(Factory_T()),
    i(BEptr->i_ref_GAMBIT()), 
    d(BEptr->d_ref_GAMBIT()), 
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
    }

T_gambit::T_gambit(Abstract_T *in) : 
    BEptr(in),
    i(BEptr->i_ref_GAMBIT()), 
    d(BEptr->d_ref_GAMBIT()),
    member_variable(false) 
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
    }

T_gambit::T_gambit(const T_gambit &in) : 
    BEptr(in.BEptr->pointerCopy_GAMBIT()),
    i(BEptr->i_ref_GAMBIT()), 
    d(BEptr->d_ref_GAMBIT()),
    member_variable(false) 
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
    }

T_gambit& T_gambit::operator=( const T_gambit& in ) 
{
    if (this != &in)
    {
        BEptr->pointerAssign_GAMBIT(in.BEptr);
    }
    return *this;
}   

T_gambit::~T_gambit()
{
    if(member_variable==false)
    {
        delete BEptr;
    }
}




//
// Class X_gambit
//

class X_gambit
{    
    private:
        bool member_variable;
    public:
        Abstract_X *BEptr;   
        T_gambit t;      

        T_gambit getT(){return T_gambit(BEptr->getT());}
        void setT(T_gambit& t_in){BEptr->setT(*t_in.BEptr);}
               
        X_gambit();
        X_gambit(Abstract_X*);
        X_gambit(const X_gambit&);
        X_gambit& operator=( const X_gambit& in );
        ~X_gambit();

        void _set_member(bool in) { member_variable = in; }
        void _print_member() { std::cout << " -- member_variable = " << member_variable << std::endl; }
};

// X_gambit class functions
X_gambit::X_gambit() : 
    BEptr(Factory_X()), 
    t(&(BEptr->t_ref_GAMBIT())),
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        t._set_member(true);        
    }

X_gambit::X_gambit(Abstract_X *in) : 
    BEptr(in), 
    t(&(BEptr->t_ref_GAMBIT())),
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        t._set_member(true);        
    }

X_gambit::X_gambit(const X_gambit &in) : 
    BEptr(in.BEptr->pointerCopy_GAMBIT()), 
    t(&(BEptr->t_ref_GAMBIT())),
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        t._set_member(true);        
    }

X_gambit& X_gambit::operator=( const X_gambit& in ) 
{
    if (this != &in)
    {
        BEptr->pointerAssign_GAMBIT(in.BEptr);
    }
    return *this;
}           
 
X_gambit::~X_gambit()
{
    if(member_variable==false)
    {
        delete BEptr;
    }
}



//
// Class Container_gambit (templated)
//

template <typename Type>
class Container_gambit {};


// Specialization: Container_gambit<int>

template <>
class Container_gambit<int>
{    
    private:
        bool member_variable;
    public:
        Abstract_Container<int> *BEptr;   
        int var;      
        void printMsg(){BEptr->printMsg();};  
               
        Container_gambit<int>();
        Container_gambit<int>(Abstract_Container<int>*);
        Container_gambit<int>(const Container_gambit<int>&);
        Container_gambit<int>& operator=( const Container_gambit<int>& in );
        ~Container_gambit<int>();

        void _set_member(bool in) { member_variable = in; }
        void _print_member() { std::cout << " -- member_variable = " << member_variable << std::endl; }
};

// Container_gambit<int> class functions

Container_gambit<int>::Container_gambit(): 
    BEptr(Factory_Container_int()),
    var(BEptr->var_ref_GAMBIT()), 
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
    }

Container_gambit<int>::Container_gambit(Abstract_Container<int> *in) : 
    BEptr(in),
    var(BEptr->var_ref_GAMBIT()), 
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
    }

Container_gambit<int>::Container_gambit(const Container_gambit<int> &in) :
    BEptr(in.BEptr->pointerCopy_GAMBIT()), 
    var(BEptr->var_ref_GAMBIT()), 
    member_variable(false) 
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
    }

Container_gambit<int>& Container_gambit<int>::operator=( const Container_gambit<int>& in ) 
{
    if (this != &in)
    {
        BEptr->pointerAssign_GAMBIT(in.BEptr);
    }
    return *this;
}   

Container_gambit<int>::~Container_gambit<int>()
{
    if(member_variable==false)
    {
        delete BEptr;
    }
}



// Specialization: Container_gambit<X_gambit>

template <>
class Container_gambit<X_gambit>
{    
    private:
        bool member_variable;
    public:
        Abstract_Container<Abstract_X> *BEptr;   
        X_gambit var;      
        void printMsg(){BEptr->printMsg();};  
               
        Container_gambit<X_gambit>();
        Container_gambit<X_gambit>(Abstract_Container<Abstract_X>*);       
        Container_gambit<X_gambit>(const Container_gambit<X_gambit>&);
        Container_gambit<X_gambit>& operator=( const Container_gambit<X_gambit>& in );
        ~Container_gambit<X_gambit>();

        void _set_member(bool in) { member_variable = in; }
        void _print_member() { std::cout << " -- member_variable = " << member_variable << std::endl; }
};

// Container_gambit<X_gambit> class functions

Container_gambit<X_gambit>::Container_gambit(): 
    BEptr(Factory_Container_X()),
    var(&(BEptr->var_ref_GAMBIT())), 
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        var._set_member(true);        
    }

Container_gambit<X_gambit>::Container_gambit(Abstract_Container<Abstract_X> *in) : 
    BEptr(in),
    var(&(BEptr->var_ref_GAMBIT())), 
    member_variable(false) 
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        var._set_member(true);        
    }

Container_gambit<X_gambit>::Container_gambit(const Container_gambit<X_gambit> &in) : 
    BEptr(in.BEptr->pointerCopy_GAMBIT()), 
    var(&(BEptr->var_ref_GAMBIT())), 
    member_variable(false) 
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        var._set_member(true);        
    }

Container_gambit<X_gambit>& Container_gambit<X_gambit>::operator=( const Container_gambit<X_gambit>& in ) 
{
    if (this != &in)
    {
        BEptr->pointerAssign_GAMBIT(in.BEptr);
    }
    return *this;
}   

Container_gambit<X_gambit>::~Container_gambit<X_gambit>()
{
    if(member_variable==false)
    {
        delete BEptr;
    }
}



// Specialization: Container_gambit<T_gambit>

template <>
class Container_gambit<T_gambit>
{    
    private:
        bool member_variable;
    public:
        Abstract_Container<Abstract_T> *BEptr;   
        T_gambit var;      
        void printMsg(){BEptr->printMsg();};  
               
        Container_gambit<T_gambit>();
        Container_gambit<T_gambit>(Abstract_Container<Abstract_T>*);
        Container_gambit<T_gambit>(const Container_gambit<T_gambit>&);
        Container_gambit<T_gambit>& operator=( const Container_gambit<T_gambit>& in );
        ~Container_gambit<T_gambit>();

        void _set_member(bool in) { member_variable = in; }
        void _print_member() { std::cout << " -- member_variable = " << member_variable << std::endl; }
};

// Container_gambit<T_gambit> class functions

Container_gambit<T_gambit>::Container_gambit(): 
    BEptr(Factory_Container_T()),
    var(&(BEptr->var_ref_GAMBIT())), 
    member_variable(false)
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        var._set_member(true);        
    }


Container_gambit<T_gambit>::Container_gambit(Abstract_Container<Abstract_T> *in) : 
    BEptr(in),
    var(&(BEptr->var_ref_GAMBIT())), 
    member_variable(false) 
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        var._set_member(true);        
    }

Container_gambit<T_gambit>::Container_gambit(const Container_gambit<T_gambit> &in) : 
    BEptr(in.BEptr->pointerCopy_GAMBIT()), 
    var(&(BEptr->var_ref_GAMBIT())), 
    member_variable(false) 
    {
        // run _set_member(true) on any member variables that are themselves wrapped classes
        var._set_member(true);        
    }

Container_gambit<T_gambit>& Container_gambit<T_gambit>::operator=( const Container_gambit<T_gambit>& in ) 
{
    if (this != &in)
    {
        BEptr->pointerAssign_GAMBIT(in.BEptr);
    }
    return *this;
}   

Container_gambit<T_gambit>::~Container_gambit<T_gambit>()
{
    if(member_variable==false)
    {
        delete BEptr;
    }
}



#endif /* __GAMBIT_WRAPPER_CLASSES_HPP__ */ 