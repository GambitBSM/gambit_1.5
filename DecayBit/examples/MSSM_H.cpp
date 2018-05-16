/**
   @file
   @brief Basic check of MSSM Higgs width
   
   @code
   g++ -std=c++11 -o MSSM_H MSSM_H.cpp -I../include/gambit/DecayBit/
   @endcode
*/

#include "MSSM_H.hpp"
#include <iostream>


int main() {
  std::array<double, 4> m{{10., 999., 999., 999.}};
  std::array<std::array<double, 4>, 4> Z;
  Z[0][0] = 1.;
  Z[0][1] = 1.;
  Z[0][2] = 1.;
  Z[0][3] = 1.;
  std::cout << MSSM_H::gamma_h_chi_0(0, 0, m, Z, 0.1) << std::endl;
}
