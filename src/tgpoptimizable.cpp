/**************************************************************************************************
 *  File:    tgpoptimizable.cpp                                                                   *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/

/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include "ublas_string.hpp"

#include "tgpoptimizable.hpp"


namespace bayesopt {

/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getBoundingBox                                                                   *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
void TGPOptimizable::getBoundingBox(vectord& lower, vectord& upper)
{
    lower = lower_bound;
    upper = upper_bound;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getUncertaintyMatrix                                                             *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
matrixd TGPOptimizable::getUncertaintyMatrix(double std_dev = 0.01)
{
    matrixd result = boost::numeric::ublas::zero_matrix<double>(dim, dim);

    for (uint element = 0; element < dim; element += 1)
    {
        result(element, element) = std_dev * std_dev;
    }

    return result;
}



/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getUncertaintyMatrixNormalized                                                   *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
matrixd TGPOptimizable::getUncertaintyMatrixNormalized(double std_dev = 0.01)
{
    matrixd result  = boost::numeric::ublas::zero_matrix<double>(dim, dim);

    for (uint element = 0; element < dim; element += 1)
    {
        result(element, element)  = std_dev / (upper_bound[element] - lower_bound[element]);
        result(element, element) *= result(element, element);
    }

    return result;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getOptParams                                                                     *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
 void TGPOptimizable::getOptParams(TgpParameters& tgp_params, Parameters& opt_params )
{
    // Initialize bopt_params
    opt_params                   = initialize_parameters_to_default();
    opt_params.n_iter_relearn    = 1;
    opt_params.init_method       = 1;
    opt_params.l_type            = L_MCMC;
    opt_params.sigma_s           = 1;
    opt_params.verbose_level     = 1;

    opt_params.kernel.name       = "kMaternARD5";

    opt_params.crit_name         = "cBEI";
    opt_params.crit_params       = vectord(2);
    opt_params.crit_params[0]    = 1.0; // Exp
    opt_params.crit_params[1]    = 0.0; // Bias

    if (opt_params.mean.coef_mean.size() == 0)
    {
        opt_params.mean.coef_mean = vectord(1,0);
    }
    else
    {
        opt_params.mean.coef_mean = vectord(opt_params.mean.coef_mean.size(), 0);
    }

    // Initialize tgp_parameters
    tgp_params.dimensions        = dim;
    tgp_params.mcmc_particles    = 10;
    tgp_params.min_data_per_leaf = 10;
    tgp_params.wheight_power     =  1;
    tgp_params.wheight_threshold =  0.00;

    setUEIMatrix(opt_params);
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setUEIMatrix                                                                     *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
void TGPOptimizable::setUEIMatrix(Parameters& opt_params)
{
    matrixd px = getUncertaintyMatrix();

    UnscentedExpectedImprovement::convertMatrixToParams(opt_params, px);
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: uniformSampling                                                                  *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
void TGPOptimizable::uniformSampling(vectord& yy, vecOfvec& xx, uint& points_per_dim)
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
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: unnormalizeVector                                                                *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
void TGPOptimizable::unnormalizeVector(vectord& x)
{
    x = boost::numeric::ublas::element_prod(upper_bound - lower_bound, x) + lower_bound;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: unnormalizeVector                                                                *
 *  Class      : TGPOptimizable                                                                   *
 **************************************************************************************************/
void TGPOptimizable::normalizeVector(vectord& x)
{
    x = boost::numeric::ublas::element_div(x - lower_bound, upper_bound - lower_bound);
}


/**
 * Converts object of TGPOptimizable to Json::Value
 *
 * @return  Json value that describes the object metadata
 */
Json::Value TGPOptimizable::getJson(void)
{
    Json::Value output;
                output["name"       ] = Json::Value(name);
                output["ymin"       ] = Json::Value(ymin);
                output["ymax"       ] = Json::Value(ymax);
                output["dim"        ] = Json::Value(dim );
                output["lower_bound"] = Json::Value(bayesopt::utils::ublas_toString(lower_bound));
                output["upper_bound"] = Json::Value(bayesopt::utils::ublas_toString(upper_bound));

    return output;
}


/**
 * Loads Json metadata into TGPOptimizable object
 *
 * @param config Metadata in form of Json::Value
 */
void TGPOptimizable::loadJson(Json::Value config)
{
    if (config["name"].isNull() != true)
        name = config["name"].asString();

    if (config["ymin"].isNull() != true)
        ymin = config["ymin"].asDouble();

    if (config["ymax"].isNull() != true)
        ymax = config["ymax"].asDouble();

    if (config["dim"].isNull() != true)
        dim = config["dim"].asUInt64();

    if (config["lower_bound"].isNull() != true)
        lower_bound = bayesopt::utils::string_toUblasVectord(config["lower_bound"].asString());

    if (config["upper_bound"].isNull() != true)
        upper_bound = bayesopt::utils::string_toUblasVectord(config["upper_bound"].asString());
}

} // End of namespace bayesopt
