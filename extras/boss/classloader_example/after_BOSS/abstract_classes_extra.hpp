// Forward declarations
template <typename T1>
class Abstract_Container;

class Abstract_T;
class Abstract_X;

class T;
class X;


// Template specializations depending on 
template <> // This specialization is only needed to go from <X> to <Abstract_X>
class Abstract_Container<X> : public Abstract_Container<Abstract_X>  
{};

template <> // This specialization is only needed to go from <T> to <Abstract_T>
class Abstract_Container<T> : public Abstract_Container<Abstract_T>
{};
