#ifndef __ABSTRACT_CLASSES_EXTRA_HPP__
#define __ABSTRACT_T_HPP__

#include "abstract_Container.hpp"

class T;
class X;

template <> // This specialization is only needed to go from <X> to <Abstract_X>
class Abstract_Container<X> : public Abstract_Container<Abstract_X>  
{};

template <> // This specialization is only needed to go from <T> to <Abstract_T>
class Abstract_Container<T> : public Abstract_Container<Abstract_T>
{};

#endif  /* __ABSTRACT_CLASSES_EXTRA_HPP__ */
