//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Variadic utilty functions.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Gregory Martinez
///          (gregory.david.martinez@gmail.com)
///  \date Aug 2014
///
///  *********************************************

#ifndef TYPE_INDEX_HPP
#define TYPE_INDEX_HPP

#include <typeinfo>
#include <functional>

namespace Gambit
{
        struct type_index
        {
        public:
                friend struct type_equal_to;
                
                type_index(const std::type_info& __rhs) noexcept
                : _M_target(&__rhs) { }

                bool operator==(const Gambit::type_index& __rhs) const noexcept
                { return *_M_target == *__rhs._M_target; }

                bool operator!=(const Gambit::type_index& __rhs) const noexcept
                { return *_M_target != *__rhs._M_target; }

                bool operator<(const Gambit::type_index& __rhs) const noexcept
                { return _M_target->before(*__rhs._M_target); }

                bool operator<=(const Gambit::type_index& __rhs) const noexcept
                { return !__rhs._M_target->before(*_M_target); }

                bool operator>(const Gambit::type_index& __rhs) const noexcept
                { return __rhs._M_target->before(*_M_target); }

                bool operator>=(const Gambit::type_index& __rhs) const noexcept
                { return !_M_target->before(*__rhs._M_target); }

                size_t hash_code() const noexcept
                { return _M_target->hash_code(); }

                const char* name() const
                { return _M_target->name(); }

        private:
                const std::type_info* _M_target;
        };
        
        struct type_hasher 
        {
                std::size_t operator()(const Gambit::type_index& code) const
                {
                        return code.hash_code();
                }
        };
        
        struct type_equal_to 
        {
                bool operator()(const Gambit::type_index& lhs, const Gambit::type_index& rhs) const
                {
                        return *lhs._M_target == *rhs._M_target;
                }
        };
}

namespace std
{
        template<typename _Tp> struct hash;
        template<typename _Tp> struct equal_to;
        
        template<>
        struct hash<Gambit::type_index>
        {
                size_t operator()(const Gambit::type_index& __ti) const noexcept
                { return __ti.hash_code(); }
        };
        
        template<>
        struct equal_to<Gambit::type_index>
        {
                size_t operator()(const Gambit::type_index& lhs, const Gambit::type_index& rhs) const noexcept
                { return lhs == rhs; }
        };
}

#endif
