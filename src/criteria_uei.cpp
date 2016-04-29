/**************************************************************************************************
 *  File:    criteria_uei.cpp                                                                     *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/

#include "criteria/criteria_uei.hpp"
#include "criteria/criteria_ei.hpp"
#include "prob_distribution.hpp"
#include "ublas_cholesky.hpp"
#include "ublas_extra.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "log.hpp"

namespace bayesopt
{

/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
UnscentedExpectedImprovement::UnscentedExpectedImprovement(void)
{
    _criteria = NULL;
}

/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
UnscentedExpectedImprovement::UnscentedExpectedImprovement(uint dim)
{
    _criteria = NULL;
    _dim      = dim;
    _scale    = 1;
    _alpha    = 0.0;
}

/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: Destructor                                                                       *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
UnscentedExpectedImprovement::~UnscentedExpectedImprovement(void)
{
    if (_criteria != NULL) delete _criteria;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: init                                                                             *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::init(NonParametricProcess* proc)
{
    mProc  = proc;
    _scale = 1;
    _alpha = 0.0;
    _dim   = proc -> getDim();
    _px    = boost::numeric::ublas::identity_matrix<double>(_dim) * 0.01 * std::sqrt(_dim + _scale);

    if (_criteria == NULL)
    {
        _criteria = new BiasedExpectedImprovement();
        _criteria -> init(proc);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setParameters                                                                    *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::setParameters(const vectord& params)
{
    if (_criteria != NULL) _criteria -> setParameters(params);

    if (params(2)  < 0) _scale = 3 - _dim;
    else                _scale = params(2);

    if (params(3) >= 0) _alpha = params(3);

    setUncertaintyMatrix(params);
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setUncertaintyMatrix                                                             *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::setUncertaintyMatrix(const vectord& params)
{
    // Uncertainty Matrix Initialization
    _px = matrixd(_dim, _dim);

    for     (uint row = 0; row < _dim; row += 1)
    {
        for (uint col = 0; col < _dim; col += 1)
        {
            _px(row, col) = params( 4 + col + (row * _dim) );
        }
    }

    // Scale Matrix
    _px *= (_dim + _scale);

    // Square root Matrix
    if (!isDiag(_px))
    {
        matrixd L;

        utils::cholesky_decompose(_px, L);

        _px = L;
    }
    else
    {
        for (uint index = 0; index < _px.size1(); index += 1)
        {
            _px(index, index) = std::sqrt(_px(index, index));
        }
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: operator()                                                                       *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::getSamples(const vectord& x, std::vector<vectord>& xx, std::vector<double>& w)
{
    xx.clear();
    w .clear();
    xx.push_back(x);
    w .push_back(_scale / (_dim + _scale));

    // Calculate query_i
    for (uint column = 0; column < _px.size2(); column += 1)
    {
        xx.push_back(x - boost::numeric::ublas::column(_px, column));
        xx.push_back(x + boost::numeric::ublas::column(_px, column));
        w .push_back(0.5 / (_dim + _scale));
        w .push_back(0.5 / (_dim + _scale));
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: operator()                                                                       *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
double UnscentedExpectedImprovement::operator() (const vectord& x)
{
    double               e_uei = 0.0;
    double               d_uei = 0.0;
    std::vector<double > ei_i;
    std::vector<vectord> xx;
                         xx.push_back(x);

    // Calculate query_i
    for (uint column = 0; column < _px.size2(); column += 1)
    {
        xx.push_back(x - boost::numeric::ublas::column(_px, column));
        xx.push_back(x + boost::numeric::ublas::column(_px, column));
    }

    // Calculate EI_i and E[UEI]
    for (uint index = 0; index < xx.size(); index += 1)
    {
        // Calculate ei_i according to underlying _criteria
        ei_i.push_back( (*_criteria)(xx[index]) );

        FILE_LOG(logDEBUG) << "cUEI Query " << index << " = " << xx[index];

        // e_uei = sum_{i=0}^{2_dim+1} (ei_i * w_i)
        if (index == 0) e_uei += ei_i[index] * ( _scale / (_dim + _scale) );
        else            e_uei += ei_i[index] * (  0.5   / (_dim + _scale) );
    }

    // Normalize e_uei
    e_uei /= xx.size();

    // Calculate V[UEI]
    for (uint index = 0; index < xx.size(); index += 1)
    {
        // d_uei = sum_{i=0}^{2_dim+1} ( (ei_i-e_uei)^2 * w_i)
        if (index == 0) d_uei += ( ei_i[index] - e_uei ) * ( ei_i[index] - e_uei ) * ( _scale / (_dim + _scale) );
        else            d_uei += ( ei_i[index] - e_uei ) * ( ei_i[index] - e_uei ) * (  0.5   / (_dim + _scale) );
    }

    // Normalize d_uei
    d_uei  = std::sqrt(d_uei / xx.size());

    FILE_LOG(logDEBUG) << "e_uei = " << e_uei;
    FILE_LOG(logDEBUG) << "d_uei = " << d_uei;

    // Return criteria
    return ( ((1 -_alpha) * e_uei) - (_alpha * d_uei) );
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: isDiag                                                                           *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
bool UnscentedExpectedImprovement::isDiag(matrixd matrix)
{
    if (matrix.size1() == matrix.size2())
    {
        for (    uint row = 0; row < matrix.size1(); row += 1)
        {
            for (uint col = 0; col < matrix.size2(); col += 1)
            {
                if (row != col)
                {
                    if (matrix(row, col) != 0.0)
                    {
                        return false;
                    }
                }
            }
        }

        return true;
    }
    else
    {
        return false;
    }

    return false;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: convertMatrixToParams                                                            *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::convertMatrixToParams(bopt_params& params, matrixd& px)
{
    if (px.size1() != px.size2()) return;

    uint dim = px.size1();

    for     (uint row = 0; row < dim; row += 1)
    {
        for (uint col = 0; col < dim; col += 1)
        {
            params.crit_params[4 + col + (row * dim)] = px(row, col);
        }
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: convertMatrixToParams                                                            *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::convertMatrixToParams(Parameters& params, matrixd& px)
{
    if (px.size1() != px.size2()) return;

    params.input.noise_matrix = px;
}

} //namespace bayesopt
