//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Type definition header for utilities for module CosmoBit.
///
///  Compile-time registration of type definitions
///  required for the rest of the code to
///  communicate with CosmoBit.
///
///  Add to this if you want to define a new type
///  for the functions in CosmoBit to return, but
///  you don't expect that type to be needed by
///  any other modules.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 Mar
///
///  *********************************************


#ifndef __CosmoBit_utils_hpp__
#define __CosmoBit_utils_hpp__

#include <valarray>
#include <vector>
#include <map>
#include <string>
#include <sstream>

namespace Gambit
{

  namespace CosmoBit
  {


    class fast_interpolation {
        // fast interpolation for grids defined on equally-spaced log space
        private:
            int grid_size;
            double Delta_logx;
            std::valarray<double> x_grid;
            std::valarray<double> y_grid;
        public:
            fast_interpolation(std::valarray<double>& x_grid0, std::valarray<double>& y_grid0)
            {
                x_grid = x_grid0;
                y_grid = y_grid0;
                grid_size = x_grid.size();
                Delta_logx = (log(x_grid[grid_size-1]) - log(x_grid[0]))/(grid_size-1);
            }
            double interp(double x)
            {
                if (x <= x_grid[0])
                {return y_grid[0];}
                if (x >= x_grid[grid_size-1])
                {return y_grid[grid_size-1];}

                double intpart_d;
                double fracpart = std::modf((log(x) - log(x_grid[0]))/Delta_logx, &intpart_d);
                int intpart = lround(intpart_d);

                return y_grid[intpart] * (1 - fracpart) + y_grid[intpart+1]*fracpart;
            }
    };

    double entropy_density_SM(double T, bool T_in_eV=false);

    std::vector<double> set_nu_masses(double mNu1, double mNu2, double mNu3, int& N_ncdm);
    std::string m_ncdm_classInput(std::map<std::string,double> NuMasses_SM);
  }
}

#endif // defined __CosmoBit_types_hpp__
