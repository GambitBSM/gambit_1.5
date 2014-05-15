#ifndef __ABSTRACT_CLASSES_HPP__
#define __ABSTRACT_CLASSES_HPP__
#include <iostream>


// MOVED TO: abstract_classes_extra.hpp
// class T;

class Abstract_T
{
    public:

        virtual void printMe() {std::cout << "Called virtual function" << std::endl;}

        virtual int& i_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_T* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
        virtual void pointerAssign_GAMBIT(Abstract_T *in){std::cout << "Called virtual function" << std::endl;}     

        virtual double& d_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
};



// MOVED TO: abstract_classes_extra.hpp
// class X;

class Abstract_X
{
    public:

        // Class methods
        virtual Abstract_T& t_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_T* getT_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
        Abstract_T* getT()
        {
            return getT_GAMBIT();
        }

        virtual void setT_GAMBIT(Abstract_T& t_in) {std::cout << "Called virtual function" << std::endl;};
        void setT(Abstract_T& t_in)
        {
            setT_GAMBIT(t_in);
        }

        virtual Abstract_X* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void pointerAssign_GAMBIT(Abstract_X *in){std::cout << "Called virtual function" << std::endl;}  
};



template <typename T1>
class Abstract_Container
{};


template <>
class Abstract_Container<Abstract_X>
{
    public:

        virtual Abstract_X& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_X>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<Abstract_X> *in){std::cout << "Called virtual function" << std::endl;}  
};

// MOVED TO: abstract_classes_extra.hpp
// template <> // This specialization is only needed to go from <X> to <Abstract_X>
// class Abstract_Container<X> : public Abstract_Container<Abstract_X>  
// {};



template <>
class Abstract_Container<int>
{
    public:

        virtual int& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<int>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<int> *in){std::cout << "Called virtual function" << std::endl;}  
};



template <> 
class Abstract_Container<Abstract_T>
{
    public:

        virtual Abstract_T& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_T>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<Abstract_T> *in){std::cout << "Called virtual function" << std::endl;}  
};

// MOVED TO: abstract_classes_extra.hpp
// template <> // This specialization is only needed to go from <T> to <Abstract_T>
// class Abstract_Container<T> : public Abstract_Container<Abstract_T>
// {};


#endif  /* __ABSTRACT_CLASSES_HPP__ */
