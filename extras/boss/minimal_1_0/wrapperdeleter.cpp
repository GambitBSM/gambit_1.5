#include "wrapperdeleter.hpp"

void wrapper_deleter(nspace3::Y__BOSS* wptr)
{
    delete wptr;
}

void wrapper_deleter(nspace1::nspace2::X__BOSS* wptr)
{
    delete wptr;
}
