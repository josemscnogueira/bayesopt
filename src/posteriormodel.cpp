/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>

   BayesOpt is free software: you can redistribute it and/or modify it
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/
#include <cfloat>

#include "log.hpp"
#include "posteriormodel.hpp"
#include "posterior_fixed.hpp"
#include "posterior_empirical.hpp"
#include "posterior_mcmc.hpp"
#include "criteria/criteria_uei.hpp"

namespace bayesopt
{

  PosteriorModel* PosteriorModel::create(size_t dim, Parameters params, randEngine& eng)
  {
    switch (params.l_type)
      {
      case L_FIXED: return new PosteriorFixed(dim,params,eng);
      case L_EMPIRICAL: return new EmpiricalBayes(dim,params,eng);
      case L_DISCRETE: // TODO:return new
      case L_MCMC: return new MCMCModel(dim,params,eng);
      case L_ERROR:
      default:
	throw std::invalid_argument("Learning type not supported");
      }
  };

  PosteriorModel::PosteriorModel(size_t dim, Parameters parameters,
				 randEngine& eng):
    mParameters(parameters), mDims(dim), mMean(dim, parameters)
  {}

  PosteriorModel::~PosteriorModel()
  { } // Default destructor


  void PosteriorModel::setSamples(const matrixd &x, const vectord &y)
  {
    mData.setSamples(x,y);
    mMean.setPoints(mData.mX);  //Because it expects a vecOfvec instead of a matrixd
  }

  void PosteriorModel::setSamples(const matrixd &x)
  {
    mData.setSamples(x);
    mMean.setPoints(mData.mX);  //Because it expects a vecOfvec instead of a matrixd
  }

  void PosteriorModel::setSamples(const vectord &y)
  {
    mData.setSamples(y);
  }


  void PosteriorModel::setSample(const vectord &x, double y)
  {
    matrixd xx(1,x.size());  vectord yy(1);
    row(xx,0) = x;           yy(0) = y;
    mData.setSamples(xx,yy);
    mMean.setPoints(mData.mX);  //Because it expects a vecOfvec instead of a matrixd
  }

  void PosteriorModel::addSample(const vectord &x, double y)
  {  mData.addSample(x,y); mMean.addNewPoint(x);  };


/**
 *
 */
double PosteriorModel::getValueUnscented(const vectord& x)
{
    // Get sigma points
    std::vector<vectord> xx;
    std::vector<double>  w;
    UnscentedExpectedImprovement::getSigmaPoints(x                                ,
                                                 mParameters.input.unscented_scale,
                                                 mDims                            ,
                                                 mParameters.input.noise_matrix   ,
                                                 xx                               ,
                                                 w                                );

    // Calculate expected observations
    std::vector<double> yy(xx.size());
    for (size_t idx_sigma = 0; idx_sigma < xx.size(); ++idx_sigma)
    {
        yy[idx_sigma] = this->getPrediction(xx[idx_sigma])->getMean();
    }

    // Calculate expected outcome
    return std::inner_product(yy.begin(), yy.end(), w.begin(), 0.0);
}


/**
 *
 */
double PosteriorModel::getValueAtMinimumUnscented(void)
{
    return getValueUnscented(getPointAtMinimumUnscented());
}


/**
 *
 */
vectord PosteriorModel::getPointAtMinimumUnscented(void)
{
    int    sample_min_idx   = -1;
    double sample_min_value = DBL_MAX;

    for (size_t sample = 0; sample < mData.getNSamples(); ++sample)
    {
        const vectord  x = mData.getSampleX(sample);
        const double  uy = getValueUnscented(x);

        if (uy < sample_min_value)
        {
            sample_min_value = uy;
            sample_min_idx   = sample;
        }
    }

    if ( (sample_min_idx >= 0                  ) &&
         (sample_min_idx <  mData.getNSamples())   )
    {
        return mData.getSampleX(sample_min_idx);
    }

    return vectord();
}


} //namespace bayesopt
