#ifndef __Objects__
#define __Objects__

#include <string>

/* Objects Of Interest (OOI) are specified by "TagType". The various bits,
//     backends, modules, and the core can then easily refer to these 
//     variables by the TagType. The type of variable can also be found 
//     easily via the TagType. ---Abram  
// NOTE: Once these macros are fully ready, neither the users nor anyone in
//     the collaboration needs to change them. Thus, to see how useful these
//     macros are, please read only the demo files. */

namespace Tags {
    template <typename TagType> 
    inline std::string tagName() {
        /* A function to return a string for the Tag name. */
        std::string errorMessage = "Totally undeclared TagType with typeid: \n  ";
        errorMessage += typeid(TagType).name();
        errorMessage += "\n";
        return errorMessage;
    }
}

template <typename TagType>
struct ObjectType {
/* a nice way to fetch what type an object has, based only on its tag.
//    for compatibility with "backend" idea, name changed to ObjectType. */
    typedef double type;  // type for a numeric observable, by default.
};

template <typename TagType>
struct ObjectOfInterest {
/* An observable, likelihood, or any other such "object of interest".
//      It is easily categorized in terms of a "TagType". */
    typedef typename ObjectType<TagType>::type ValueType;
    static ValueType *value;
};


#define ASSOCIATE_TAG_OBJECT(TAG_TYPE, VALUE_TYPE) \
/* Create a specification of ObjectType, and immediately
//     associate the TAG_TYPE with the ObjectOfInterest's VALUE_TYPE.
//     Also, make a function which returns a nice string for the TAG_TYPE. */ \
namespace Tags { \
    struct TAG_TYPE{}; \
    template<> \
    inline std::string tagName<Tags::TAG_TYPE>() { \
        return #TAG_TYPE; \
    } \
} \
\
template<> \
struct ObjectType<Tags::TAG_TYPE> { \
    typedef VALUE_TYPE type; \
}; \
\
template<> \
ObjectOfInterest<Tags::TAG_TYPE>::ValueType \
                * ObjectOfInterest<Tags::TAG_TYPE>::value = NULL;


#define PASS_ADDRESS(TAG_TYPE, ADDRESS) \
/* Pass the address of an actual variable in memory to the ObjectOfInterest. */ \
ObjectOfInterest<Tags::TAG_TYPE>::value = ADDRESS; 


#define OOI(TAG_TYPE) \
/* Conveniencee macro for the user to get the value of the ObjectOfInterest. */ \
*ObjectOfInterest<Tags::TAG_TYPE>::value


#endif /* defined(__Objects__) */
