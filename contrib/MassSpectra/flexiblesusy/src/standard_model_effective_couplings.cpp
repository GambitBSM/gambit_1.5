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

#include "standard_model_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "standard_model.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {
namespace standard_model {

using namespace effective_couplings;

#define INPUTPARAMETER(parameter) model.get_input().parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define DERIVEDPARAMETER(parameter) model.##parameter()
#define PHASE(parameter) model.get_##parameter()
#define PHYSICAL(parameter) model.get_physical().parameter

Standard_model_effective_couplings::Standard_model_effective_couplings(
   const Standard_model& model_,
   const softsusy::QedQcd& qedqcd_,
   const Physical_input& input_)
   : model(model_), qedqcd(qedqcd_), physical_input(input_)
   , rg_improve(true), include_qcd_corrections(true)
   , Vd(PHYSICAL(Vd)), Ud(PHYSICAL(Ud)), Vu(PHYSICAL(Vu)), Uu(PHYSICAL(Uu)), Ve
      (PHYSICAL(Ve)), Ue(PHYSICAL(Ue)), ZZ(PHYSICAL(ZZ))

   , eff_CphhVPVP(0), eff_CphhVGVG(0), eff_CpAhVPVP(0), eff_CpAhVGVG(0)

{
}

void Standard_model_effective_couplings::calculate_effective_couplings()
{
   const standard_model::Standard_model sm(initialise_SM());
   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_parameters(model.get());

   const double saved_mt = PHYSICAL(MFu(2));
   PHYSICAL(MFu(2)) = qedqcd.displayPoleMt();

   const auto Mhh = PHYSICAL(Mhh);
   run_SM_strong_coupling_to(sm, 0.5 * Mhh);
   calculate_eff_CphhVPVP();
   run_SM_strong_coupling_to(sm, Mhh);
   calculate_eff_CphhVGVG();

   PHYSICAL(MFu(2)) = saved_mt;
   model.set_scale(scale);
   model.set(saved_parameters);

}

void Standard_model_effective_couplings::set_model(const Standard_model& model_)
{
   model = model_;
   copy_mixing_matrices_from_model();
}

void Standard_model_effective_couplings::copy_mixing_matrices_from_model()
{
   Vd = PHYSICAL(Vd);
   Ud = PHYSICAL(Ud);
   Vu = PHYSICAL(Vu);
   Uu = PHYSICAL(Uu);
   Ve = PHYSICAL(Ve);
   Ue = PHYSICAL(Ue);
   ZZ = PHYSICAL(ZZ);

}

standard_model::Standard_model Standard_model_effective_couplings::initialise_SM() const
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_input);

   sm.initialise_from_input(qedqcd);

   return sm;
}

void Standard_model_effective_couplings::run_SM_strong_coupling_to(standard_model::Standard_model sm, double m)
{
   sm.run_to(m);
   model.set_g3(sm.get_g3());
}

std::complex<double> Standard_model_effective_couplings::scalar_scalar_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      if (m_loop > m_decay) {
         result = 1 + 0.06754745576155852*Sqr(g3);
      }

   }

   return result;
}

std::complex<double> Standard_model_effective_couplings::scalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         scalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

std::complex<double> Standard_model_effective_couplings::pseudoscalar_fermion_qcd_factor(double m_decay, double m_loop) const
{
   std::complex<double> result(1.0, 0.0);

   if (include_qcd_corrections) {
      const auto g3 = MODELPARAMETER(g3);
      result = 1.0 + 0.025330295910584444*Sqr(g3) *
         pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);

   }

   return result;
}

double Standard_model_effective_couplings::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double Standard_model_effective_couplings::scalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(23.75 - 1.1666666666666667*Nf)*
      Sqr(g3);
   const double nnlo_qcd = 0.000641623890917771*Power(g3,4)*(370.1956513893174
      + 2.375*l + (-47.18640261449638 + 0.6666666666666666*l)*Nf +
      0.9017702481178881*Sqr(Nf));
   const double nnnlo_qcd = 0.000016252523020247696*Power(g3,6)*(467.683620788
      + 122.440972222*l + 10.9409722222*Sqr(l));

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double Standard_model_effective_couplings::pseudoscalar_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));

   const auto g3 = MODELPARAMETER(g3);

   const double nlo_qcd = 0.025330295910584444*(24.25 - 1.1666666666666667*Nf)*
      Sqr(g3);
   const double nnlo_qcd = 0.000641623890917771*Power(g3,4)*(171.54400563089382
      + 5*l);
   const double nnnlo_qcd = 0;

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

double Standard_model_effective_couplings::get_hhVPVP_partial_width() const
{
   const double mass = PHYSICAL(Mhh);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CphhVPVP);
}

double Standard_model_effective_couplings::get_hhVGVG_partial_width() const
{
   const double mass = PHYSICAL(Mhh);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CphhVGVG);
}

double Standard_model_effective_couplings::get_AhVPVP_partial_width() const
{
   const double mass = PHYSICAL(MAh);
   return 0.0049735919716217296 * Power(mass, 3.0) * AbsSqr(eff_CpAhVPVP);
}

double Standard_model_effective_couplings::get_AhVGVG_partial_width() const
{
   const double mass = PHYSICAL(MAh);
   return 0.039788735772973836 * Power(mass, 3.0) * AbsSqr(eff_CpAhVGVG);
}

std::complex<double> Standard_model_effective_couplings::CpFdhhbarFdPL(int gt1, int gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_467;
   std::complex<double> tmp_468;
   for (int j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_469;
      std::complex<double> tmp_470;
      for (int j1 = 0; j1 < 3; ++j1) {
         tmp_470 += Conj(Ud(gt3,j1))*Yd(j1,j2);
      }
      tmp_469 += tmp_470;
      tmp_468 += (Conj(Vd(gt1,j2))) * tmp_469;
   }
   tmp_467 += tmp_468;
   result += (-0.7071067811865475) * tmp_467;

   return result;
}

std::complex<double> Standard_model_effective_couplings::CpFuhhbarFuPL(int gt1, int gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_471;
   std::complex<double> tmp_472;
   for (int j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_473;
      std::complex<double> tmp_474;
      for (int j1 = 0; j1 < 3; ++j1) {
         tmp_474 += Conj(Uu(gt3,j1))*Yu(j1,j2);
      }
      tmp_473 += tmp_474;
      tmp_472 += (Conj(Vu(gt1,j2))) * tmp_473;
   }
   tmp_471 += tmp_472;
   result += (0.7071067811865475) * tmp_471;

   return result;
}

std::complex<double> Standard_model_effective_couplings::CpFehhbarFePL(int gt1, int gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_475;
   std::complex<double> tmp_476;
   for (int j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_477;
      std::complex<double> tmp_478;
      for (int j1 = 0; j1 < 3; ++j1) {
         tmp_478 += Conj(Ue(gt3,j1))*Ye(j1,j2);
      }
      tmp_477 += tmp_478;
      tmp_476 += (Conj(Ve(gt1,j2))) * tmp_477;
   }
   tmp_475 += tmp_476;
   result += (-0.7071067811865475) * tmp_475;

   return result;
}

double Standard_model_effective_couplings::CphhVWpconjVWp() const
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   double result = 0.0;

   result = 0.5*v*Sqr(g2);

   return result;
}

std::complex<double> Standard_model_effective_couplings::CpAhFdbarFdPL(int gt2, int gt3) const
{
   const auto Yd = MODELPARAMETER(Yd);

   std::complex<double> result;

   std::complex<double> tmp_479;
   std::complex<double> tmp_480;
   for (int j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_481;
      std::complex<double> tmp_482;
      for (int j1 = 0; j1 < 3; ++j1) {
         tmp_482 += Conj(Ud(gt3,j1))*Yd(j1,j2);
      }
      tmp_481 += tmp_482;
      tmp_480 += (Conj(Vd(gt2,j2))) * tmp_481;
   }
   tmp_479 += tmp_480;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_479;

   return result;
}

std::complex<double> Standard_model_effective_couplings::CpAhFubarFuPL(int gt2, int gt3) const
{
   const auto Yu = MODELPARAMETER(Yu);

   std::complex<double> result;

   std::complex<double> tmp_483;
   std::complex<double> tmp_484;
   for (int j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_485;
      std::complex<double> tmp_486;
      for (int j1 = 0; j1 < 3; ++j1) {
         tmp_486 += Conj(Uu(gt3,j1))*Yu(j1,j2);
      }
      tmp_485 += tmp_486;
      tmp_484 += (Conj(Vu(gt2,j2))) * tmp_485;
   }
   tmp_483 += tmp_484;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_483;

   return result;
}

std::complex<double> Standard_model_effective_couplings::CpAhFebarFePL(int gt2, int gt3) const
{
   const auto Ye = MODELPARAMETER(Ye);

   std::complex<double> result;

   std::complex<double> tmp_487;
   std::complex<double> tmp_488;
   for (int j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_489;
      std::complex<double> tmp_490;
      for (int j1 = 0; j1 < 3; ++j1) {
         tmp_490 += Conj(Ue(gt3,j1))*Ye(j1,j2);
      }
      tmp_489 += tmp_490;
      tmp_488 += (Conj(Ve(gt2,j2))) * tmp_489;
   }
   tmp_487 += tmp_488;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_487;

   return result;
}

void Standard_model_effective_couplings::calculate_eff_CphhVPVP()
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto MVWp = PHYSICAL(MVWp);
   const auto decay_mass = PHYSICAL(Mhh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * scalar_fermion_qcd_factor(decay_mass,
         MFd(gI1)) * CpFdhhbarFdPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * scalar_fermion_qcd_factor(decay_mass,
         MFu(gI1)) * CpFuhhbarFuPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFehhbarFePL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   result += -0.5 * CphhVWpconjVWp() * vev * AS1(decay_scale / Sqr(MVWp)) / Sqr
      (MVWp);


   result *= 0.18926819071273507 * physical_input.get(
      Physical_input::alpha_em_0) * Sqrt(qedqcd.displayFermiConstant());

   eff_CphhVPVP = result;

}

void Standard_model_effective_couplings::calculate_eff_CphhVGVG()
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(Mhh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFdhhbarFdPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpFuhhbarFuPL(gI1, gI1) * vev * AS12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= std::complex<double>(0.75,0.);

   if (include_qcd_corrections) {
      result *= scalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant()
      );

   eff_CphhVGVG = result;

}

void Standard_model_effective_couplings::calculate_eff_CpAhVPVP()
{
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const auto MFe = PHYSICAL(MFe);
   const auto decay_mass = PHYSICAL(MAh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 0.3333333333333333 * pseudoscalar_fermion_qcd_factor(
         decay_mass, MFd(gI1)) * CpAhFdbarFdPL(gI1, gI1) * vev * AP12(decay_scale
         / Sqr(MFd(gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += 1.3333333333333333 * pseudoscalar_fermion_qcd_factor(
         decay_mass, MFu(gI1)) * CpAhFubarFuPL(gI1, gI1) * vev * AP12(decay_scale
         / Sqr(MFu(gI1))) / MFu(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFebarFePL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFe(
         gI1))) / MFe(gI1);
   }
   result *= std::complex<double>(2.0,0.);


   result *= 0.18926819071273507 * physical_input.get(
      Physical_input::alpha_em_0) * Sqrt(qedqcd.displayFermiConstant());

   eff_CpAhVPVP = result;

}

void Standard_model_effective_couplings::calculate_eff_CpAhVGVG()
{
   const auto g3 = MODELPARAMETER(g3);
   const auto MFd = PHYSICAL(MFd);
   const auto MFu = PHYSICAL(MFu);
   const double alpha_s = 0.07957747154594767*Sqr(g3);
   const auto decay_mass = PHYSICAL(MAh);
   const auto decay_scale = 0.25 * Sqr(decay_mass);

   const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   std::complex<double> result = 0;
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFdbarFdPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFd(
         gI1))) / MFd(gI1);
   }
   for (int gI1 = 0; gI1 < 3; ++gI1) {
      result += CpAhFubarFuPL(gI1, gI1) * vev * AP12(decay_scale / Sqr(MFu(
         gI1))) / MFu(gI1);
   }
   result *= std::complex<double>(1.5,0.);

   if (include_qcd_corrections) {
      result *= pseudoscalar_scaling_factor(decay_mass);
   }


   result *= 0.12617879380849006 * alpha_s * Sqrt(qedqcd.displayFermiConstant()
      );

   eff_CpAhVGVG = result;

}

} // namespace standard_model
} // namespace flexiblesusy
