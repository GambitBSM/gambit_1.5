// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef PARALLEL_H
#define PARALLEL_H

#include "config.h"

#ifdef ENABLE_THREADS

#include <future>
#include <utility>

namespace flexiblesusy {

template<typename F, typename... Ts>
inline std::future<typename std::result_of<F(Ts...)>::type>
run_async(F&& f, Ts&&... params)
{
   return std::async(std::launch::async,
                     std::forward<F>(f),
                     std::forward<Ts>(params)...);
}

} // namespace flexiblesusy

#endif

#endif
