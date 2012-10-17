#include <boost/math/distributions/students_t.hpp> // for student t distribution

#include "student_t_process.hpp"
#include "cholesky.hpp"
#include "trace_ublas.hpp"

  
StudentTProcess::StudentTProcess(double noise):
  NonParametricProcess(noise)
{}  // Constructor



StudentTProcess::~StudentTProcess()
{} // Default destructor




double StudentTProcess::negativeLogLikelihood(size_t index)
{
  matrixd K = computeCorrMatrix();
  size_t n = K.size1();
  
  matrixd L(n,n);
  cholesky_decompose(K,L);

  vectord colU = (*mMean)(mGPXX);;

  //TODO: Replace by transform
  //  for (size_t ii=0; ii< n; ii++) 
  //  colU(ii) = 

  vectord alphU(colU);
  boost::numeric::ublas::inplace_solve(L,alphU,boost::numeric::ublas::lower_tag());
  double eta = inner_prod(colU,alphU);
  
  vectord alphY(mGPY);
  boost::numeric::ublas::inplace_solve(L,alphY,boost::numeric::ublas::lower_tag());
  double mu     = inner_prod(colU,alphY) / eta;
  double YInvRY = inner_prod(mGPY,alphY);
    
  double sigma = (YInvRY - mu*mu*eta) / (n-1);

  double negloglik = 0.5*( (n-1)*log(sigma) + trace(L) + log(eta) );

  return negloglik;
}


int StudentTProcess::prediction( const vectord &query,
				 double& yPred, double& sPred)
{
  size_t n = mGPXX.size();
  vectord rInvR(n);
  double kn;
  double uInvRr, rInvRr;
  double meanf = mMean->getMean(query);
  
  vectord colR = computeCrossCorrelation(query);
  kn = (*mKernel)(query, query);
  
  noalias(rInvR) = prod(colR,mInvR);	
  rInvRr = inner_prod(rInvR,colR);
  uInvRr = inner_prod(mUInvR,colR);
  
  svectord colMu(n,mMu);
  vectord yumu = mGPY - meanf*colMu;
  
  yPred = meanf*mMu + inner_prod( rInvR, yumu );
  sPred = sqrt( mSig * (kn - rInvRr + (meanf - uInvRr) * (meanf - uInvRr) 
			/ mUInvRUDelta ) );

  return n-1;
}
	

int StudentTProcess::precomputePrediction()
{
  size_t nSamples = mGPXX.size();
  vectord colU = (*mMean)(mGPXX);

  //TODO: Replace by transform
  //  for (size_t ii=0; ii< nSamples; ii++) 
  //  colU(ii) = (*mMean)(mGPXX[ii]);

  mUInvR = prod(colU,mInvR);
  mUInvRUDelta = inner_prod(mUInvR,colU);
  
  vectord YInvR(nSamples);
  double YInvRY;
  
  mMu =  inner_prod(mUInvR,mGPY) / mUInvRUDelta;
  
  noalias(YInvR) = prod(mGPY,mInvR);
  YInvRY = inner_prod(YInvR,mGPY);
  
  mSig = (YInvRY - mMu*mMu*mUInvRUDelta) / (nSamples-1);

  return 1;
}


	
double StudentTProcess::negativeExpectedImprovement(const vectord& query,
						    size_t g)
{
  size_t dof = mGPXX.size() - 1;
  boost::math::students_t d(dof);

  double yPred,sPred;
  double yMin = getValueAtMinimum();
  prediction(query,yPred,sPred);

  double yDiff = yMin - yPred; 
  double yNorm = yDiff / sPred;
  
  if (g != 1)
    {
      std::cout << "Students t EI with exponent not yet supported." << std::endl;
      return 0.0;
    }
  else
    {
      return -( yDiff * cdf(d,yNorm) + 
		(dof*sPred+yNorm*yDiff)/(dof-1) * pdf(d,yNorm) );
    }
  
}  // negativeExpectedImprovement

double StudentTProcess::lowerConfidenceBound(const vectord& query,
					     double beta)
{    
  size_t n = mGPXX.size();
  double yPred,sPred;
  prediction(query,yPred,sPred);
  return yPred - beta*sPred/sqrt(n);
}  // lowerConfidenceBound

double StudentTProcess::negativeProbabilityOfImprovement(const vectord& query,
							 double epsilon)
{  
  size_t dof = mGPXX.size() - 1;
  boost::math::students_t d(dof);
  double yPred,sPred;
  double yMin = getValueAtMinimum();
  prediction(query,yPred,sPred);

  return -cdf(d,(yMin - yPred + epsilon)/sPred);
}  // negativeProbabilityOfImprovement