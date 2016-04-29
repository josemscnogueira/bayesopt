#ifndef __TGPPARAMETERS_H__
#define __TGPPARAMETERS_H__

#include <vector>
#include <string>


typedef struct
{
    uint                     dimensions;
    uint                     mcmc_particles;
    uint                     min_data_per_leaf;
    double                   wheight_power;
    double                   wheight_threshold;
    double                   samples_to_save;
    std::vector<std::string> others;
}
tgp_parameters;


#endif // __TGPPARAMETERS_H__
