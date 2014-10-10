#include "wrapperdeleter.hpp"

void wrapper_deleter(X_GAMBIT* wptr)
{
    delete wptr;
}

void wrapper_deleter(Y_GAMBIT* wptr)
{
    delete wptr;
}
