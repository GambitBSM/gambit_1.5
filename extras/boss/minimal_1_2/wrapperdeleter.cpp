#include "wrapperdeleter.hpp"

void wrapper_deleter(X__BOSS* wptr)
{
    delete wptr;
}

void wrapper_deleter(Y__BOSS* wptr)
{
    delete wptr;
}
