#ifndef __NULLPTR_CHECK__
#define __NULLPTR_CHECK__

#include <stdexcept>

template<typename T>
T nullptr_check(T ptr)
{
	if(!ptr)
	{
		throw std::runtime_error("BOSS says: You are trying to dereference an uninitialized pointer. We shall have null of it.");
	}

	return ptr;
}

#endif /* __NULLPTR_CHECK__ */