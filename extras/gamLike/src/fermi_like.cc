// GamLike interface to *bayes
//
// Christoph Weniger, 2014-07-13
// <c.weniger@uva.nl>

#include "gamLike.h"
#include <gsl/gsl_integration.h>
#include <iostream>

namespace gamLike
{
    // Global pointer to dwarfs and gc object
    extern GamCombLike * dwarfsCombLike;
    extern GamCombLike * gcCombLike;

    extern "C" double fermi_dwarfs_likelihood_(double & mass, double & scaling, double & sv, double * Energy, double * dNdE)
    {
        std::vector<double> E;
        std::vector<double> dPhidE;
        for ( int i = 0; i < 100; i++)
        {
            E.push_back(Energy[i]);
            dPhidE.push_back(dNdE[i] * scaling * scaling * sv / mass / mass / 8 / 3.141593);
        }
        return dwarfsCombLike->lnL(E, dPhidE);
    }

    extern "C" double fermi_dwarfs_likelihood_vec_(const std::vector<double> & E, const std::vector<double> & dPhidE)
    {
        return dwarfsCombLike->lnL(E, dPhidE);
    }

    extern "C" double fermi_gc_likelihood_(double & mass, double & scaling, double & sv, double * Energy, double * dNdE)
    {
        std::vector<double> E;
        std::vector<double> dPhidE;
        for ( int i = 0; i < 100; i++)
        {
            E.push_back(Energy[i]);
            dPhidE.push_back(dNdE[i] * scaling * scaling * sv / mass / mass / 8 / 3.141593);
        }
        return gcCombLike->lnL(E, dPhidE);
    }

    extern "C" double fermi_gc_likelihood_vec_(const std::vector<double> & E, const std::vector<double> & dPhidE)
    {
        return gcCombLike->lnL(E, dPhidE);
    }
}
