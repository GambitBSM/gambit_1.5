#ifndef __ABSTRACT_CLASSES_EXTRA_HPP__
#define __ABSTRACT_T_HPP__
// Forward declarations:

// MOVED FROM: abstract_T.hpp
class T;
// CLASS DEFINITION IN: abstract_T.hpp
class Abstract_T;

// MOVED FROM: abstract_X.hpp
class X;
// CLASS DEFINITION IN: abstract_X.hpp
class Abstract_X;

// CLASS DEFINITION IN: abstract_Container.hpp
template <typename T1>
class Abstract_Container;

// Specializations which require the above forward declarations:

// MOVED FROM: abstract_Container.hpp
template <>
class Abstract_Container<Abstract_X>
{
    public:

        virtual Abstract_X& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_X>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<Abstract_X> *in){std::cout << "Called virtual function" << std::endl;}  
};

// MOVED FROM: abstract_Container.hpp
template <> 
class Abstract_Container<Abstract_T>
{
    public:

        virtual Abstract_T& var_ref_GAMBIT() {std::cout << "Called virtual function" << std::endl;};

        virtual void printMsg() {std::cout << "Called virtual function" << std::endl;};

        virtual Abstract_Container<Abstract_T>* pointerCopy_GAMBIT() {std::cout << "Called virtual function" << std::endl;};
 
        virtual void pointerAssign_GAMBIT(Abstract_Container<Abstract_T> *in){std::cout << "Called virtual function" << std::endl;}  
};

// MOVED FROM: abstract_Container.hpp
template <> // This specialization is only needed to go from <X> to <Abstract_X>
class Abstract_Container<X> : public Abstract_Container<Abstract_X>  
{};
// MOVED FROM: abstract_Container.hpp
template <> // This specialization is only needed to go from <T> to <Abstract_T>
class Abstract_Container<T> : public Abstract_Container<Abstract_T>
{};

#endif  /* __ABSTRACT_CLASSES_EXTRA_HPP__ */
