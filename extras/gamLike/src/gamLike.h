//
//           *** GamLike ***
//
// Gamma-ray likelihoods from tables
//
// Author: Christoph Weniger <c.weniger@uva.nl>
// Date: 2014-07-13
// Version: 0.1
//

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <gsl/gsl_integration.h>

namespace gamLike
{
    typedef std::vector<std::pair<std::vector<double>, std::vector<double> > > LnL_vs_flux;
    typedef std::vector<std::pair<double, double> > E_bins;


    //////////////////////////////////////////
    // Marginalized likelihood for one target
    //////////////////////////////////////////

    class GamLike: public gsl_function
    {
        public:
            GamLike(std::string filename, double J_mean, double J_err, int method, int filetype);
            double lnL(std::vector<double> Phi);
            E_bins getEbins();

        private:
            static double invoke(double J, void* par);
            int loadData_lnLtab(std::string filename);
            int loadData_Sigma(std::string filename);
            double lnL(std::vector<double> Phi, double J);
            double lnL_lnLtab(std::vector<double> Phi, double J);
            double lnL_Sigma(std::vector<double> Phi, double J);
            int method;  // 0: mean value; 1: marginalization; 2: profiling; 3: profiling with 1/J factor
            int filetype;  // 0: tabulated lnL; 1: flux with covariance matrix

            double J_mean;  // log10(GeV^2 cm^-5 sr)
            double J_err;   // log10(GeV^2 cm^-5 sr)
            gsl_integration_workspace * gsl_workspace;  // GSL workspace
            std::vector<double> my_Phi, flux;
            std::vector<std::vector<double> > Sigma;
            std::string filename;  // filename
            E_bins e_bins;  // Energy bins [GeV]
            LnL_vs_flux lnL_vs_flux;  // list of lnL as function of flux in each bin [GeV/cm^2/s]
    };


    ////////////////////////////////////////////////
    // Marginalized likelihood for multiple targets
    ////////////////////////////////////////////////

    class GamCombLike
    {
        public:
            GamCombLike(int method): method(method) {};
            void add_target(std::string filename, double J_mean, double J_err, int filetype);
            double lnL(std::vector<double> Phi);
            double lnL(std::vector<double> E, std::vector<double> dPhidE);
            E_bins getEbins();

        private:
            std::vector<GamLike> my_targets;
            E_bins e_bins;
            int method;
    };


    ////////////////////
    // Helper functions 
    ////////////////////

    // Normal distribution with mean x0 and standard deviation dx
    double normal_pdf(double x0, double dx, double x);

    // Interpolation in log-log space
    double interpolate(double x, const std::vector<double> & xlist, const std::vector<double> & ylist, bool zerobound);

    class IntSpec: public gsl_function
    {
        public:
            IntSpec(double* Energy, double* dNdE);
            IntSpec(std::vector<double> Energy, std::vector<double> dNdE);
            static double invoke(double x, void* params);
            double integrate(double x0, double x1);

        private:
            std::vector<double> Energy;
            std::vector<double> dNdE;
            gsl_integration_workspace * gsl_workspace;  // GSL workspace
    };
}
