/**************************************************************************************************
 *  File:    tgpoptimizable.hpp                                                                   *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/
#ifndef __TGPOPTIMIZABLE_HPP__
#define __TGPOPTIMIZABLE_HPP__


/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include <string>

#include "specialtypes.hpp"
#include "tgpparameters.h"
#include "parameters.h"
#include "criteria/criteria_uei.hpp"


/**************************************************************************************************
 *  Interface: TGPOptimizable                                                                     *
 **************************************************************************************************/
class TGPOptimizable
{
public:
    // Attributes
    std::string name;
    double      ymax;
    double      ymin;
    uint        dim;
    vectord     lower_bound;
    vectord     upper_bound;

    // Destructor
    virtual ~TGPOptimizable(void) {};

    // Methods
    virtual double evaluate      (vectord   x                   ) = 0;
    virtual void   getBoundingBox(vectord& lower, vectord& upper)
    {
        lower = lower_bound;
        upper = upper_bound;
    };

    virtual matrixd getUncertaintyMatrix(double std_dev = 0.01)
    {
        matrixd result = boost::numeric::ublas::zero_matrix<double>(dim, dim);

        for (uint element = 0; element < dim; element += 1)
        {
            result(element, element) = std_dev * std_dev;
        }

        return result;
    };

    virtual matrixd getUncertaintyMatrixNormalized(double std_dev = 0.01)
    {
        matrixd result  = boost::numeric::ublas::zero_matrix<double>(dim, dim);

        for (uint element = 0; element < dim; element += 1)
        {
            result(element, element)  = std_dev / (upper_bound[element] - lower_bound[element]);
            result(element, element) *= result(element, element);
        }

        return result;
    }


    virtual void getOptParams(tgp_parameters& tgp_params, bopt_params& opt_params )
    {
        // Initialzie bopt_params
        opt_params                   = initialize_parameters_to_default();
        opt_params.n_iter_relearn    = 1;
        opt_params.init_method       = 1;
        opt_params.l_type            = L_MCMC;
        opt_params.sigma_s           = 1;
        opt_params.verbose_level     = 1;
        opt_params.kernel.name       = "kMaternARD5";
        opt_params.crit_name         = "cBEI";
        opt_params.crit_params[0]    = 1;   // exp
        opt_params.crit_params[1]    = 0.00; // bias
        opt_params.n_crit_params     = 2;
        for (uint index = 0; index < opt_params.mean.n_coef; index += 1)
        {
            opt_params.mean.coef_mean[index] = 0;
        }
        // Initialzie tgp_parameters
        tgp_params.dimensions        = dim;
        tgp_params.mcmc_particles    = 10;
        tgp_params.min_data_per_leaf = 10;
        tgp_params.wheight_power     =  1;
        tgp_params.wheight_threshold =  0.00;
        tgp_params.samples_to_save   =  5;
    };

    virtual void setUEIMatrix(bopt_params& opt_params)
    {
        matrixd px = getUncertaintyMatrix();

        bayesopt::UnscentedExpectedImprovement::convertMatrixToParams(opt_params, px);
    };

    virtual void uniformSampling(vectord& yy, vecOfvec& xx, uint& points_per_dim)
    {
        while (std::pow(points_per_dim, dim) > 100000) points_per_dim /= 2;

        xx.clear();
        yy = vectord(std::pow(points_per_dim, dim), 0.0);

        for (uint sample = 0; sample < yy.size(); sample += 1)
        {
            vectord query = vectord(dim, 0.0);

            for (uint index = 0; index < dim; index += 1)
            {
                uint x_i  = sample % ( (uint) std::pow(points_per_dim, index + 1) );
                     x_i /=          ( (uint) std::pow(points_per_dim, index + 0) );

                query(index) = lower_bound(index) + ( (1.0 / points_per_dim) * x_i * (upper_bound(index) - lower_bound(index)));
            }

            yy(sample) = evaluate(query);
            xx.push_back         (query);
        }
    };

    virtual void unnormalizeVector(vectord& x)
    {
        for (uint index = 0; index < dim; index += 1)
        {
            x[index] = ((upper_bound[index] - lower_bound[index]) * x[index]) + lower_bound[index];
        }
    }

    virtual void normalizeVector(vectord& x)
    {
        for (uint index = 0; index < dim; index += 1)
        {
            x[index] = (x[index] - lower_bound[index]) / (upper_bound[index] - lower_bound[index]);
        }
    }
};


#endif // _TGPOPTIMIZABLE_H_
