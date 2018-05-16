/**
   @file
   @brief Basic check of MSSM Higgs width

   @code
   g++ -std=c++11 -o MSSM_H MSSM_H.cpp -I../include/gambit/DecayBit/
   @endcode

   The example matches an SLHA file with the blocks

   @code
   BLOCK MASS  # Mass Spectrum
   25     125.             # h
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

   BLOCK ALPHA
   -8.46891780E-01   # Mixing angle in the neutral Higgs boson sector
   @endcode

   for which I find this result from HDECAY:

   @code
   DECAY        25     1.02087441E+05   # h decays
   5.43844340E-03    2     1000022   1000022   # BR(h -> ~chi_10 ~chi_10)
   5.04017313E-02    2     1000023   1000023   # BR(h -> ~chi_20 ~chi_20)
   2.01559758E-05    2     1000025   1000025   # BR(h -> ~chi_30 ~chi_30)
   4.44948403E-02    2     1000035   1000035   # BR(h -> ~chi_40 ~chi_40)
   3.42370645E-02    2     1000022   1000023   # BR(h -> ~chi_10 ~chi_20)
   5.00810176E-05    2     1000022   1000025   # BR(h -> ~chi_10 ~chi_30)
   7.21398155E-02    2     1000022   1000035   # BR(h -> ~chi_10 ~chi_40)
   2.36593112E-05    2     1000023   1000025   # BR(h -> ~chi_20 ~chi_30)
   1.39380325E-01    2     1000023   1000035   # BR(h -> ~chi_20 ~chi_40)
   3.92722704E-03    2     1000025   1000035   # BR(h -> ~chi_30 ~chi_40)
   @endcode

*/

#include "MSSM_H.hpp"
#include <iostream>


int main() {
  std::array<double, 4> m{{10., 20., 30., 40.}};
  std::array<std::array<double, 4>, 4> Z;
  Z[0][0] = 9.86364430E-01;
  Z[0][1] = -5.31103553E-02;
  Z[0][2] = 1.46433995E-01;
  Z[0][3] = -5.31186117E-02;
  Z[1][0] =  9.93505358E-02;
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
      std::cout << MSSM_H::gamma_h_chi_0(i, j, m, Z, -8.46891780E-01)
                << " " << i + 1
                << " " << j + 1 <<  std::endl;
    }
  }
}
