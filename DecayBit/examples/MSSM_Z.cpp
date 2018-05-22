/**
   @file
   @brief Basic check of MSSM Z width

   @code
   g++ -std=c++11 -o MSSM_Z MSSM_Z.cpp -I../include/gambit/DecayBit/
   @endcode

   The example matches an SLHA file with the blocks

   @code
   BLOCK MASS  # Mass Spectrum
   1000022     10.              # ~chi_10
   1000023     20.              # ~chi_20
   1000025     30.              # ~chi_30
   1000035     40.              # ~chi_40

   BLOCK NMIX  # Neutralino Mixing Matrix
   1  1     9.86364430E-01   # N_11
   1  2    -5.31103553E-02   # N_12
   1  3     1.46433995E-01   # N_13
   1  4    -5.31186117E-02   # N_14
   2  1     9.93505358E-02   # N_21
   2  2     9.44949299E-01   # N_22
   2  3    -2.69846720E-01   # N_23
   2  4     1.56150698E-01   # N_24
   3  1    -6.03388002E-02   # N_31
   3  2     8.77004854E-02   # N_32
   3  3     6.95877493E-01   # N_33
   3  4     7.10226984E-01   # N_34
   4  1    -1.16507132E-01   # N_41
   4  2     3.10739017E-01   # N_42
   4  3     6.49225960E-01   # N_43
   4  4    -6.84377823E-01   # N_44
   @endcode
   
   for which I find this result from microMEGAs:

   @code
   Z :   total width=2.51E+00 
   and Branchings:
   6.701431E-02 Z -> ne,Ne
   3.380425E-02 Z -> e,E
   6.701431E-02 Z -> nm,Nm
   3.380425E-02 Z -> m,M
   6.701431E-02 Z -> nl,Nl
   3.372794E-02 Z -> l,L
   1.162932E-01 Z -> u,U
   1.495033E-01 Z -> d,D
   1.495033E-01 Z -> s,S
   1.162464E-01 Z -> c,C
   1.487689E-01 Z -> b,B
   2.388469E-05 Z -> ~o1,~o1
   1.125856E-04 Z -> ~o1,~o2
   1.057151E-04 Z -> ~o2,~o2
   1.779237E-03 Z -> ~o1,~o3
   5.688950E-03 Z -> ~o2,~o3
   1.367656E-05 Z -> ~o3,~o3
   2.241813E-04 Z -> ~o1,~o4
   1.723165E-04 Z -> ~o2,~o4
   9.195154E-03 Z -> ~o3,~o4
   @endcode
*/

#include "MSSM_Z.hpp"
#include <iostream>

int main() {
  std::array<double, 4> m_0{{10., 20., 30., 40.}};
  std::array<std::array<double, 4>, 4> Z;
  Z[0][0] = 9.86364430E-01;
  Z[0][1] = -5.31103553E-02;
  Z[0][2] = 1.46433995E-01;
  Z[0][3] = -5.31186117E-02;
  Z[1][0] = 9.93505358E-02;
  Z[1][1] = 9.44949299E-01;
  Z[1][2] = -2.69846720E-01;
  Z[1][3] = 1.56150698E-01;
  Z[2][0] = -6.03388002E-02;
  Z[2][1] = 8.77004854E-02;
  Z[2][2] = 6.95877493E-01;
  Z[2][3] = 7.10226984E-01;
  Z[3][0] = -1.16507132E-01;
  Z[3][1] = 3.10739017E-01;
  Z[3][2] = 6.49225960E-01;
  Z[3][3] = -6.84377823E-01;

  for (int i = 0; i <= 3; i += 1) {
    for (int j = 0; j <= 3; j += 1) {
      std::cout << MSSM_Z::gamma_chi_0(i, j, m_0, Z)
                << " " << i + 1
                << " " << j + 1 <<  std::endl;
    }
  }
}
