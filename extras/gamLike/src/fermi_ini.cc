// GamLike interface to *bayes
//
// Christoph Weniger, 2014-07-13
// <c.weniger@uva.nl>

#include "gamLike.h"
#include <iostream>

// Method:
// 0: mean value
// 1: marginalization
// 2: profiling
// 3: profiling with additional 1/J

namespace gamLike
{
    // Global pointer to dwarfs object
    GamCombLike * dwarfsCombLike;
    GamCombLike * gcCombLike;
    std::string path;

    extern "C" void set_data_path_(const char * str)
    {
        path = str;
    }

    extern "C" void fermi_dwarfs_init_(int method)
    {
        std::cout << "Initializing Fermi LAT dwarfs routines." << std::endl;

#ifdef GAMLIKE_DATA_PATH
        set_data_path_(std::string(GAMLIKE_DATA_PATH).c_str());
#endif

        dwarfsCombLike = new GamCombLike(method);
        dwarfsCombLike->add_target(path+"like_bootes_I.txt", 18.8, 0.22, 0);
        dwarfsCombLike->add_target(path+"like_canes_venatici_I.txt", 17.7, 0.26, 0);
        dwarfsCombLike->add_target(path+"like_canes_venatici_II.txt", 17.9, 0.25, 0);
        dwarfsCombLike->add_target(path+"like_carina.txt", 18.1, 0.23, 0);
        dwarfsCombLike->add_target(path+"like_coma_berenices.txt", 19.0, 0.25, 0);
        dwarfsCombLike->add_target(path+"like_draco.txt", 18.8, 0.16, 0);
        dwarfsCombLike->add_target(path+"like_fornax.txt", 18.2, 0.21, 0);
        dwarfsCombLike->add_target(path+"like_hercules.txt", 18.1, 0.25, 0);
        dwarfsCombLike->add_target(path+"like_leo_I.txt", 17.7, 0.18, 0);
        dwarfsCombLike->add_target(path+"like_leo_II.txt", 17.6, 0.18, 0);
        dwarfsCombLike->add_target(path+"like_leo_IV.txt", 17.9, 0.28, 0);
        dwarfsCombLike->add_target(path+"like_sculptor.txt", 18.6, 0.18, 0);
        dwarfsCombLike->add_target(path+"like_segue_1.txt", 19.5, 0.29, 0);
        dwarfsCombLike->add_target(path+"like_sextans.txt", 18.4, 0.27, 0);
        dwarfsCombLike->add_target(path+"like_ursa_major_I.txt", 18.3, 0.24, 0);
        dwarfsCombLike->add_target(path+"like_ursa_major_II.txt", 19.3, 0.28, 0);
        dwarfsCombLike->add_target(path+"like_ursa_minor.txt", 18.8, 0.19, 0);
        dwarfsCombLike->add_target(path+"like_willman_1.txt", 19.1, 0.31, 0);
    }

    extern "C" void fermi_gc_init_(int method, int type)
    {
        std::cout << "Initializing Fermi GC routines." << std::endl;

#ifdef GAMLIKE_DATA_PATH
        set_data_path_(std::string(GAMLIKE_DATA_PATH).c_str());
#endif

        gcCombLike = new GamCombLike(method);
        if (type == 0)
        {
            gcCombLike->add_target(path+"like_GC_Daylan.txt", 23.8026, 0.01, 0);
        }
        if (type == 1)
        {
            gcCombLike->add_target(path+"like_GC_Calore.txt", 23.315, 0.01, 1);
        }
    }
}
