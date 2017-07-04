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
        _criteria->init(proc);
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

    if (params(2) <  0) _scale = 3 - _dim;
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

    for (size_t row = 0; row < _dim; ++row)
    {
        for (size_t col = 0; col < _dim; ++col)
        {
            _px(row, col) = params(4 + col + (row * _dim));
        }
    }

    _px = UnscentedExpectedImprovement::convertMatrixNoise(_px, _scale, _dim);
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: operator()                                                                       *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::getSamples(const vectord& x, std::vector<vectord>& xx, std::vector<double>& w)
{
    UnscentedExpectedImprovement::getSigmaPoints(x, _scale, _dim, _px, xx, w, false);
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
    for (size_t column = 0; column < _px.size2(); ++column)
    {
        xx.push_back(x - boost::numeric::ublas::column(_px, column));
        xx.push_back(x + boost::numeric::ublas::column(_px, column));
    }

    // Calculate EI_i and E[UEI]
    for (size_t idx = 0; idx < xx.size(); ++idx)
    {
        // Calculate ei_i according to underlying _criteria
        ei_i.push_back( (*_criteria)(xx[idx]) );

        FILE_LOG(logDEBUG) << "cUEI Query " << idx << " = " << xx[idx];

        // e_uei = sum_{i=0}^{2_dim+1} (ei_i * w_i)
        if (idx == 0) e_uei += ei_i[idx] * ( _scale / (_dim + _scale) );
        else          e_uei += ei_i[idx] * (  0.5   / (_dim + _scale) );
    }

    // Normalize e_uei
    e_uei /= xx.size();

    // Calculate V[UEI]
    for (size_t idx = 0; idx < xx.size(); idx += 1)
    {
        // d_uei = sum_{i=0}^{2_dim+1} ( (ei_i-e_uei)^2 * w_i)
        if (idx == 0) d_uei += ( ei_i[idx] - e_uei ) * ( ei_i[idx] - e_uei ) * ( _scale / (_dim + _scale) );
        else          d_uei += ( ei_i[idx] - e_uei ) * ( ei_i[idx] - e_uei ) * (  0.5   / (_dim + _scale) );
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
        for (size_t row = 0; row < matrix.size1(); ++row)
        {
            for (size_t col = 0; col < matrix.size2(); ++col)
            {
                if (row != col)
                {
                    if (std::abs(matrix(row, col)) > std::numeric_limits<double>::epsilon())
                    {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    return false;
}


/**
 *
 */
matrixd UnscentedExpectedImprovement::convertMatrixNoise(const matrixd& matrix,
                                                         const double   scale ,
                                                         const int      dim   )
{
    matrixd matrix_output = matrix;

    // Scale Matrix
    matrix_output *= (dim + scale);

    // Square root Matrix
    if (!isDiag(matrix_output))
    {
        matrixd L;
        utils::cholesky_decompose(matrix_output, L);

        matrix_output = L;
    }
    else
    {
        for (size_t idx = 0; idx < matrix_output.size1(); ++idx)
        {
            matrix_output(idx, idx) = std::sqrt(matrix_output(idx, idx));
        }
    }

    return matrix_output;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: convertMatrixToParams                                                            *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::convertMatrixToParams(bopt_params& params, const matrixd px)
{
    if (px.size1() != px.size2()) return;

    size_t dim = px.size1();

    for (size_t row = 0; row < dim; ++row)
    {
        for (size_t col = 0; col < dim; ++col)
        {
            params.crit_params[4 + col + (row * dim)] = px(row, col);
            params.input.noise[0 + col + (row * dim)] = px(row, col);
        }
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: convertMatrixToParams                                                            *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::convertMatrixToParams(Parameters& params, const matrixd px)
{
    if (px.size1() != px.size2()) return;

    size_t dim = px.size1();

    params.input.noise_matrix = px;
    params.crit_params.resize(4+dim*dim, true);

    for (size_t row = 0; row < dim; ++row)
    {
        for (size_t col = 0; col < dim; ++col)
        {
            params.crit_params[4 + col + (row * dim)] = px(row, col);
        }
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getSigmaPoints                                                                   *
 *  Class      : UnscentedExpectedImprovement                                                     *
 **************************************************************************************************/
void UnscentedExpectedImprovement::getSigmaPoints(const vectord&          x             ,
                                                  const double            scale         ,
                                                  const int               dim           ,
                                                  const matrixd&          matrix_noise  ,
                                                  std::vector<vectord>&   xx            ,
                                                  std::vector<double>&    w             ,
                                                  const bool              matrix_convert)
{
    const size_t n = dim;

    assert(matrix_noise.size1() == n);
    assert(matrix_noise.size2() == n);
    assert(x.size()             == n);

    matrixd px;
    if (matrix_convert) px = UnscentedExpectedImprovement::convertMatrixNoise(matrix_noise, scale, dim);
    else                px = matrix_noise;

    // Output variable intialization
    xx = std::vector<vectord>();
    w  = std::vector<double>();
    xx.push_back(x);
    w .push_back(scale / (dim + scale));

    // Calculate query_i
    for (size_t col = 0; col < n; col += 1)
    {
        xx.push_back(x - boost::numeric::ublas::column(px, col));
        xx.push_back(x + boost::numeric::ublas::column(px, col));
        w .push_back(0.5 / (dim + scale));
        w .push_back(0.5 / (dim + scale));
    }
}

} //namespace bayesopt
