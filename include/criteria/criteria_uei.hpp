/**************************************************************************************************
 *  File:    criteria_uei.hpp                                                                     *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/

#ifndef  _CRITERIA_UEI_HPP_
#define  _CRITERIA_UEI_HPP_

#include "criteria_functors.hpp"

namespace bayesopt
{

class UnscentedExpectedImprovement: public Criteria
{
public:
    // Constructor
    UnscentedExpectedImprovement(void);
    UnscentedExpectedImprovement(uint dim);

    // Destructor
    virtual ~UnscentedExpectedImprovement(void);

    // Methods
    void        init                (NonParametricProcess* proc  );
    void        setParameters       (const vectord&        params);
    size_t      nParameters         (void                        );
    double      operator()          (const vectord&        x     );
    std::string name                (void                        )
        {return "cUEI"; };
    void        getSamples          (const vectord& x, std::vector<vectord>& xx, std::vector<double>& w);

    // Static mathods
    static void convertMatrixToParams(bopt_params& params, matrixd& px);
    static void convertMatrixToParams(Parameters&  params, matrixd& px);

private:
    // Attributes
    Criteria* _criteria;
    double    _scale;     // parameters(2)
    double    _dim;
    matrixd   _px;
    double    _alpha;     // parameters(3)

    // Methods
    bool isDiag              (const matrixd  matrix);
    void setUncertaintyMatrix(const vectord& params);
};

inline size_t UnscentedExpectedImprovement::nParameters() { return (2 + _criteria -> nParameters() + (_dim * _dim)); };


} //namespace bayesopt


#endif
