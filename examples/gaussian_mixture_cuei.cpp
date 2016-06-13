/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>

   Also, author of this file: Jose Nogueira <josemscnogueira@gmail.com>

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

#include <boost/numeric/ublas/matrix.hpp>

#include "testfunctions.hpp"
#include "param_loader.hpp"
#include "criteria/criteria_uei.hpp"

int main(int nargs, char *args[])
{
    bayesopt::Parameters par;

    par = initialize_parameters_to_default();

    // Optimization and fuction params
    par.n_init_samples    =  30;
    par.n_iterations      = 150 - par.n_init_samples;
    par.random_seed       =   0;
    par.verbose_level     =   1;
    par.noise             = 1e-10;
    par.sigma_s           =   0.0642;

    // Learning params
    par.n_iter_relearn    = 1;
    par.init_method       = 1;
    par.l_type            = L_MCMC;

    // Kernel Params
    par.kernel.name       = "kMaternARD5";

    // Learning criterion Params
    par.crit_name         = "cUEI";
    par.crit_params       = vectord(4);
    par.crit_params[0]    = 1.0;  // Exp
    par.crit_params[1]    = 0.0;  // Bias
    par.crit_params[2]    = 1.00; // scale
    par.crit_params[3]    = 0.00; // alpha

    // Create input space noise matrix
    matrixd px = boost::numeric::ublas::identity_matrix<double>(2) * 0.02; // Input space is normalized

    // Update BayesOpt params according to that matrix
    bayesopt::UnscentedExpectedImprovement::convertMatrixToParams(par, px);

    // Initialize optimization variables
    GaussianMixture2DNormalized gm_opt(par);
    vectord                     result(2);

    // Start Optimization
    gm_opt.optimize(result);

    // Print Result
    std::cout << "Result: " << result << "->"
	          << gm_opt.evaluateSample(result) << std::endl;


    return 0;
}
