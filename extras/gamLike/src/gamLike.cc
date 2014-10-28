//
//           *** GamLike ***
//
// Gamma-ray likelihoods from tables
//
// Author: Christoph Weniger <c.weniger@uva.nl>
// Date: 2014-07-13
// Version: 0.1
//


#include "gamLike.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <gsl/gsl_integration.h>

namespace gamLike
{
    double normal_pdf(double x0, double dx, double x)
    {
        return exp(-pow(x-x0, 2)/2/pow(dx, 2))/sqrt(2*3.141592)/dx;
    }

    double interpolate(double x, const std::vector<double> & xlist, const std::vector<double> & ylist, bool zerobound)
    {
        double x0, x1, y0, y1;
        int i = 0;
        if (zerobound)
        {
            if (x<xlist.front()) return 0;
            if (x>xlist.back()) return 0;
        }
        else
        {
            if (x<xlist.front()) return ylist.front();
            if (x>xlist.back()) return ylist.back();
        }
        for (; xlist[i] < x; i++) {};  // Find index
        x0 = xlist[i-1];
        x1 = xlist[i];
        y0 = ylist[i-1];
        y1 = ylist[i]; 
        return y0 + (y1-y0) * log(x/x0) / log(x1/x0);  // lin-vs-log interpolation for lnL vs flux
    }

    GamLike::GamLike(std::string filename, double J_mean, double J_err, int method, int filetype): J_mean(J_mean), J_err(J_err), filename(filename), method(method), filetype(filetype)
    {
        if (filetype == 0)
            loadData_lnLtab(filename);
        if (filetype == 1)
            loadData_Sigma(filename);
    }

    // static member function, takes void pointer as class reference
    double GamLike::invoke(double J, void* par)
    {
        GamLike* gamLike = static_cast<GamLike*>(par);
        double a1 = normal_pdf(gamLike->J_mean, gamLike->J_err, J);
        double a2 = exp(gamLike->lnL(gamLike->my_Phi, J));
        return a1*a2;
    }

    double GamLike::lnL(std::vector<double> Phi)
    {
        double epsabs = 0.;
        double epsrel = 0.001;
        size_t limit = 10000;
        double result, error;

        // Skip marginalization/profiling if likelihood low or error small
        if (this->J_err < J_mean*0.001 or lnL(Phi, J_mean) < -50 or method == 0)
        {
            return lnL(Phi, this->J_mean);
        }
        else
        {
            function = &GamLike::invoke;
            params = this;
            my_Phi = Phi;

            if (method == 1)
            {
                // Marginalization
                gsl_workspace = gsl_integration_workspace_alloc (100000);
                gsl_integration_qags(this, J_mean-5*J_err, J_mean+5*J_err, epsabs, epsrel, limit, gsl_workspace, &result, &error);
                gsl_integration_workspace_free(gsl_workspace);

                return log(result);
            }
            if (method == 2)
            {
                // Profiling
                double L, L_max = -100000;
                for (double J = J_mean - 5*J_err; J < J_mean +5*J_err; J = J + 0.01*J_err)
                {
                    L = this->invoke(J, this);
                    if (L > L_max) L_max = L;
                }
                return log(L_max);
            }
            if (method == 3)
            {
                // Profiling Fermi style
                double L, L_max = -100000;
                for (double J = J_mean - 5*J_err; J < J_mean +5*J_err; J = J + 0.01*J_err)
                {
                    // Note: Dividing by J is what is done in
                    // Ackermann:2013yva. This introduce some hidden prior
                    // dependence, but reproduces their results.
                    L = this->invoke(J, this)/pow(10, J);
                    if (L > L_max) L_max = L;
                }
                return log(L_max);
            }
        }
    }

    E_bins GamLike::getEbins()
    {
        return e_bins;
    }

    int GamLike::loadData_lnLtab(std::string filename)
    {
        double phi, lnL;
        int nlines = 0;
        std::pair<double, double> ebin, ebin_last;
        std::vector<double> phi_list, lnL_list;

        std::cout << "Load data for " + filename << "." << std::endl;
        std::ifstream in(filename.c_str(), std::ios::binary);
        if (in.fail())
        {
            std::cout << "ERROR: Failed loading " + filename + "." << std::endl;
            exit(-1);
        }
        std::string line;
        while(std::getline(in, line)) 
        {
            if (line[0] == '#') continue;
            std::stringstream ss(line);
            ss >> ebin.first;
            ebin.first /= 1e3;  // scale to GeV
            ss >> ebin.second;
            ebin.second /= 1e3;  // scale to GeV
            ss >> phi;
            phi /= 1e3;  // scale to GeV
            ss >> lnL;
            if ( nlines == 0 ) ebin_last = ebin;
            if ( ebin_last != ebin )
            {
                e_bins.push_back(ebin_last);
                ebin_last = ebin;
                lnL_vs_flux.push_back(std::pair<std::vector<double>, std::vector<double> >(phi_list, lnL_list));
                phi_list.clear();
                lnL_list.clear();
            }
            phi += 1e-20;  // Non-zero zero flux everywhere
            phi_list.push_back(phi);
            lnL_list.push_back(lnL);
            nlines += 1;
        } 
        e_bins.push_back(ebin);
        lnL_vs_flux.push_back(std::pair<std::vector<double>, std::vector<double> >(phi_list, lnL_list));
        in.close();
    }

    int GamLike::loadData_Sigma(std::string filename)
    {
        int nlines = 0;
        double tmp;
        std::pair<double, double> ebin;

        std::cout << "Load data for " + filename << "." << std::endl;
        std::ifstream in(filename.c_str(), std::ios::binary);
        if (in.fail())
        {
            std::cout << "ERROR: Failed loading " + filename + "." << std::endl;
            exit(-1);
        }
        std::string line;
        while(std::getline(in, line)) 
        {
            if (line[0] == '#') continue;
            std::stringstream ss(line);

            // Read energy bins
            ss >> ebin.first;
            ss >> ebin.second;
            e_bins.push_back(ebin);

            // Read flux
            ss >> tmp;
            this->flux.push_back(tmp);

            // Read Sigma 
            std::vector<double> tmpVec;
            while(ss >> tmp)
            {
                tmpVec.push_back(tmp);
            }
            this->Sigma.push_back(tmpVec);
        } 
        std::cout << "Flux vector length: " << this->flux.size() << std::endl;
        std::cout << "Sigma dim: " << this->Sigma.size() << " x " << this->Sigma.begin()->size() << std::endl;
        in.close();
    }

    double GamLike::lnL(std::vector<double> Phi, double J)
    {
        double lnL = 0;
        if (filetype == 0)
            lnL = lnL_lnLtab(Phi, J);
        if (filetype == 1)
            lnL = lnL_Sigma(Phi, J);
        return lnL;
    }

    double GamLike::lnL_Sigma(std::vector<double> Phi, double J)
    {
        // Phi: int(dNdE) * sv / m**2 / 8 pi [cm^3/s/GeV**2]
        // J: J-value [GeV**2/cm^5]
        // returns \Delta \ln L
        double lnL = 0;
        int N = this->flux.size();
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                lnL -= 0.5*(Phi[i] * pow(10, J) - flux[i]) * (Phi[j] * pow(10, J) - flux[j]) * Sigma[i][j];
            }
        return lnL;
    }

    double GamLike::lnL_lnLtab(std::vector<double> Phi, double J)
    {
        // Phi: int(dNdE) * sv / m**2 / 8 pi [cm^3/s/GeV**2]
        // J: J-value [GeV**2/cm^5]
        // returns \Delta \ln L
        double lnL = 0, eMean, energy_flux;
        std::vector<double>::iterator it_Phi = Phi.begin();
        E_bins::iterator it_E_bin = e_bins.begin();
        for (LnL_vs_flux::iterator it = lnL_vs_flux.begin(); it != lnL_vs_flux.end(); it++)
        {
            eMean = sqrt(it_E_bin->first * it_E_bin->second);  // [GeV]
            energy_flux = *it_Phi * pow(10, J) * eMean;  // GeV/cm^2/s
            lnL += interpolate(energy_flux, it->first, it->second, false);
            it_Phi++;
            it_E_bin++;
        }
        return lnL;
    }

    void GamCombLike::add_target(std::string filename, double J_mean, double J_err, int filetype)
    {
        if (method == 0) J_err = 0;
        GamLike d(filename, J_mean, J_err, method, filetype);
        if ( e_bins.empty() )
        {
            my_targets.push_back(d);
            e_bins = d.getEbins();
        }
        else
        {
            if ( e_bins == d.getEbins() )
            {
                my_targets.push_back(d);
            }
            else
            {
                std::cout << "ERROR: Failed loading " + filename + "." << std::endl;
                std::cout << "Incompatible energy bins." << std::endl;
                exit(-1);
            }
        }
    }

    double GamCombLike::lnL(std::vector<double> Phi)
    {
        double lnL = 0;
        for (std::vector<GamLike>::iterator it = my_targets.begin(); it != my_targets.end(); it++)
        {
            lnL += it->lnL(Phi);
        }
        return lnL;
    }

    double GamCombLike::lnL(std::vector<double> E, std::vector<double> dPhidE)
    {
        std::vector<double> Phi;
        IntSpec intSpec(E, dPhidE);
        for (E_bins::iterator it = e_bins.begin(); it != e_bins.end(); it++)
        {
            double phi = intSpec.integrate(it->first, it->second);
            Phi.push_back(phi);
        }
        return lnL(Phi);
    }

    E_bins GamCombLike::getEbins()
    {
        return e_bins;
    }

    IntSpec::IntSpec(double* Energy, double* dNdE)
    {
        function = &IntSpec::invoke;
        params=this;
        std::copy(&Energy[0], &Energy[100], back_inserter(this->Energy));
        std::copy(&dNdE[0], &dNdE[100], back_inserter(this->dNdE));
    }

    IntSpec::IntSpec(std::vector<double> Energy, std::vector<double> dNdE)
    {
        function = &IntSpec::invoke;
        params=this;
        this->Energy = Energy;
        this->dNdE = dNdE;
    }

    double IntSpec::integrate(double x0, double x1)
    {
        //return this->invoke(sqrt(x0*x1), this)*(x1-x0);
        double epsabs = 0.;
        double epsrel = 0.001;
        size_t limit = 10000;
        double result, error;
        gsl_workspace = gsl_integration_workspace_alloc (100000);
        gsl_integration_qags(this, x0, x1, epsabs, epsrel, limit, gsl_workspace, &result, &error);
        gsl_integration_workspace_free(gsl_workspace);
        return result;
    }

    double IntSpec::invoke(double x, void* params)
    {
        IntSpec* p= static_cast<IntSpec*>(params);
        double t = interpolate(x, p->Energy, p->dNdE, true);
        return t;

    }
}
