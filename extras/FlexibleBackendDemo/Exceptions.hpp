#ifndef __Exceptions__
#define __Exceptions__

/* These classes are simple exceptions which are thrown by the other
//     various classes. ---Abram  
// NOTE: Once these macros are fully ready, neither the users nor anyone in
//     the collaboration needs to change them. Thus, to see how useful these
//     macros are, please read only the demo files. */

#include <exception>
#include <string>
#include "Objects.hpp"

using namespace std;

template <typename TagType>
struct NoTagSupportException : public std::exception {
    /* a simple exception class to handle when Tags are not supported at all. */
    const char * what () const throw ()
    {
        string errorMessage = "Unsupported TagType: \n  ";
        errorMessage += Tags::tagName<TagType>();
        errorMessage += "\n";
        return errorMessage.c_str();
    }
};
 
template <typename TagType>
struct NotProvidedException : public std::exception {
    /* a simple exception class to handle when Tags are not provided by the Bit. */
    const char * what () const throw ()
    {
        string errorMessage = "That OOI is REQUIRED, not provided. Tag: \n";
        errorMessage += Tags::tagName<TagType>();
        errorMessage += "\n";
        return errorMessage.c_str();
    }
};
 
template <typename TagType>
struct NotRequiredException : public std::exception {
    /* a simple exception class to handle when Tags are not required by the Bit. */
    const char * what () const throw ()
    {
        string errorMessage = "That OOI is PROVIDED, not required. Tag: \n";
        errorMessage += Tags::tagName<TagType>();
        errorMessage += "\n";
        return errorMessage.c_str();
    }
};
 
struct NotReadyException : public std::exception {
    /* a simple exception class to handle when the Bit is not ready to calculation. */
    const char * what () const throw ()
    {
        return "Not ready!!! Provide all of the requirements.";
    }
};
 

#endif /* defined(__Exceptions__) */
