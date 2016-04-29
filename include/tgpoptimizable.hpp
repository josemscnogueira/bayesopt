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
#include "bayesopt/parameters.hpp"
#include "criteria/criteria_uei.hpp"


namespace bayesopt {

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
    virtual double  evaluate                      (vectord  x) = 0;

    virtual void    getBoundingBox                (vectord& lower,
                                                   vectord& upper);

    virtual matrixd getUncertaintyMatrix          (double  std_dev  );
    virtual matrixd getUncertaintyMatrixNormalized(double  std_dev  );

    virtual void    getOptParams                  (TgpParameters& tgp_params,
                                                   Parameters&    opt_params);
    virtual void    setUEIMatrix                  (Parameters&    opt_params);

    virtual void    uniformSampling               (vectord& yy, vecOfvec& xx,
                                                   uint&    points_per_dim);

    virtual void    unnormalizeVector             (vectord& x);
    virtual void    normalizeVector               (vectord& x);
};

} // End of namespace bayesopt


#endif // _TGPOPTIMIZABLE_H_
