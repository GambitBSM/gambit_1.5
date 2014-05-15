#ifndef __ABSTRACT_CLASSES_HPP__
#define __ABSTRACT_CLASSES_HPP__
#include <iostream>


// class T;
class Abstract_T
{
    public:

        virtual void printMe() {std::cout << "Called virtual function" << std::endl;}

        virtual int& iRef() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_T* ptr_copy() {std::cout << "Called virtual function" << std::endl;};
        virtual void ptr_assign(Abstract_T *in){std::cout << "Called virtual function" << std::endl;}     

        virtual double& dRef() {std::cout << "Called virtual function" << std::endl;};
};



// class X;
class Abstract_X
{
    public:

        // Class methods
        virtual Abstract_T& tRef() {std::cout << "Called virtual function" << std::endl;};
        // virtual Abstract_T& tRef_wrapper() {std::cout << "Called virtual function" << std::endl;};
        // Abstract_T& tRef()
        // {
        //     return tRef_wrapper();
        // }

        virtual Abstract_T* getT_wrapper() {std::cout << "Called virtual function" << std::endl;};
        Abstract_T* getT()
        {
            return getT_wrapper();
        }

        virtual void setT_wrapper(Abstract_T& t_in) {std::cout << "Called virtual function" << std::endl;};
        void setT(Abstract_T& t_in)
        {
            setT_wrapper(t_in);
        }

        virtual void setT2_wrapper(Abstract_T** t_in, Abstract_T** t_in2) {std::cout << "Called virtual function" << std::endl;};
        void setT2(Abstract_T** t_in, Abstract_T** t_in2)
        {
            setT2_wrapper(t_in, t_in2);
        }

        virtual Abstract_X* ptr_copy() {std::cout << "Called virtual function" << std::endl;};

        virtual void ptr_assign(Abstract_X *in){std::cout << "Called virtual function" << std::endl;}  
};



template <typename T1>
class Abstract_Container
{};


template <>
class Abstract_Container<Abstract_X>
{
    public:

        virtual Abstract_X& varRef() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_X>* ptr_copy() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void ptr_assign(Abstract_Container<Abstract_X> *in){std::cout << "Called virtual function" << std::endl;}  
};

// template <> // This specialization is only needed to go from <X> to <Abstract_X>
// class Abstract_Container<X> : public Abstract_Container<Abstract_X>  
// {};



template <>
class Abstract_Container<int>
{
    public:

        virtual int& varRef() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<int>* ptr_copy() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void ptr_assign(Abstract_Container<int> *in){std::cout << "Called virtual function" << std::endl;}  
};



template <> 
class Abstract_Container<Abstract_T>
{
    public:

        virtual Abstract_T& varRef() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_T>* ptr_copy() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void ptr_assign(Abstract_Container<Abstract_T> *in){std::cout << "Called virtual function" << std::endl;}  
};

// template <> // This specialization is only needed to go from <T> to <Abstract_T>
// class Abstract_Container<T> : public Abstract_Container<Abstract_T>
// {};


#endif  /* __ABSTRACT_CLASSES_HPP__ */
