#include "wrapperdeleter.hpp"

void wrapper_deleter(nspace3::Y_GAMBIT* wptr)
{
    delete wptr;
}

void wrapper_deleter(nspace1::nspace2::X_GAMBIT* wptr)
{
    delete wptr;
}
