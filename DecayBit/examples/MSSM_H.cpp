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
   1000024     10.              # ~chi_1+
   1000037     20.              # ~chi_2+

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

   BLOCK UMIX  # Chargino Mixing Matrix U
   1  1    -6.87450932E-01   # U_11
   1  2     7.26230829E-01   # U_12
   2  1     7.26230829E-01   # U_21
   2  2     6.87450932E-01   # U_22

   BLOCK VMIX  # Chargino Mixing Matrix V
   1  1    -7.26230829E-01   # V_11
   1  2     6.87450932E-01   # V_12
   2  1     6.87450932E-01   # V_21
   2  2     7.26230829E-01   # V_22

   BLOCK ALPHA
   -8.46891780E-01   # Mixing angle in the neutral Higgs boson sector
   @endcode

   for which I find this result from HDECAY:

   @code
   DECAY        25     1.02087441E+05   # h decays
   2.23325818E-01    2     1000024  -1000024   # BR(h -> ~chi_1+ ~chi_1-)
   1.94802015E-01    2     1000037  -1000037   # BR(h -> ~chi_2+ ~chi_2-)
   8.38088166E-04    2     1000024  -1000037   # BR(h -> ~chi_1+ ~chi_2-)
   8.38088166E-04    2     1000037  -1000024   # BR(h -> ~chi_2+ ~chi_1-)
   3.15536853E-03    2     1000022   1000022   # BR(h -> ~chi_10 ~chi_10)
   2.92429331E-02    2     1000023   1000023   # BR(h -> ~chi_20 ~chi_20)
   1.16944366E-05    2     1000025   1000025   # BR(h -> ~chi_30 ~chi_30)
   2.58157727E-02    2     1000035   1000035   # BR(h -> ~chi_40 ~chi_40)
   1.98642420E-02    2     1000022   1000023   # BR(h -> ~chi_10 ~chi_20)
   2.90568560E-05    2     1000022   1000025   # BR(h -> ~chi_10 ~chi_30)
   4.18553043E-02    2     1000022   1000035   # BR(h -> ~chi_10 ~chi_40)
   1.37270613E-05    2     1000023   1000025   # BR(h -> ~chi_20 ~chi_30)
   8.08680459E-02    2     1000023   1000035   # BR(h -> ~chi_20 ~chi_40)
   2.27856534E-03    2     1000025   1000035   # BR(h -> ~chi_30 ~chi_40)
   @endcode
*/

#include "MSSM_H.hpp"
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
      std::cout << MSSM_H::gamma_h_chi_0(i, j, m_0, Z, -8.46891780E-01)
                << " " << i + 1
                << " " << j + 1 <<  std::endl;
    }
  }

  std::array<double, 2> m_pm{{10., 20.}};
  std::array<std::array<double, 2>, 2> U;
  std::array<std::array<double, 2>, 2> V;
  U[0][0] = -6.87450932E-01;
  U[0][1] = 7.26230829E-01;
  U[1][0] = 7.26230829E-01;
  U[1][1] = 6.87450932E-01;

  V[0][0] = -7.26230829E-01;
  V[0][1] = 6.87450932E-01;
  V[1][0] = 6.87450932E-01;
  V[1][1] = 7.26230829E-01;

  for (int i = 0; i <= 1; i += 1) {
    for (int j = 0; j <= 1; j += 1) {
      std::cout << MSSM_H::gamma_h_chi_pm(i, j, m_pm, U, V, -8.46891780E-01)
                << " " << i + 1
                << " " << j + 1 <<  std::endl;
    }
  }
}
