/**************************************************************************************************
 *  File:    treedgaussianprocess.hpp                                                             *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/
#ifndef __TREEDGAUSSIANPROCESS_HPP__
#define __TREEDGAUSSIANPROCESS_HPP__

/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include <vector>                 // Global
#include <string>
#include <fstream>

#include "posteriormodel.hpp"     // Bayesopt
#include "prob_distribution.hpp"
#include "tgpnode.hpp"


/**************************************************************************************************
 *  Namespace: bayesopt                                                                           *
 **************************************************************************************************/
namespace bayesopt
{

/**************************************************************************************************
 *  Class: TreedGaussianProcess                                                                   *
 **************************************************************************************************/
class TreedGaussianProcess : public PosteriorModel
{
public:
    // Constructor
    TreedGaussianProcess(TgpParameters& tgp_params, Parameters& params, randEngine& eng);

    // Destructor
    ~TreedGaussianProcess();

    // Methods (redefined)
    void                     setSamples                (const matrixd &x, const vectord &y);
    void                     setSample                 (const vectord &x, double         y);
    void                     addSample                 (const vectord &x, double         y);

    void                     updateHyperParameters     (void);
    void                     fitSurrogateModel         (void);
    void                     updateSurrogateModel      (void);

    double                   evaluateCriteria          (const vectord& query);
    void                     updateCriteria            (const vectord& query);
    bool                     criteriaRequiresComparison(void);
    void                     setFirstCriterium         (void);
    bool                     setNextCriterium          (const vectord& prevResult);
    std::string              getBestCriteria           (      vectord& best);
    ProbabilityDistribution* getPrediction             (const vectord& query);

    // Methods (new)
    void                     print                     (void);
    std::vector<double>      getSplits                 (void);
    uint                     getNumberOfLeafs          (void);
    std::vector<TGPNode*>    getLeafsPtrs              (void);

private:
    // References
    TgpParameters&     tgpparams;
    randEngine&        engine;

    // Attributes
    TGPNode*           _root;
    uint               _population;
};


/**************************************************************************************************
 *  Inline methods : TreedGaussianProcess                                                         *
 **************************************************************************************************/
inline void                     TreedGaussianProcess::addSample                 (const vectord& x, double y) {        _root -> addSample                 (x, y      ); }
inline void                     TreedGaussianProcess::fitSurrogateModel         (void                      ) {        _root -> fitSurrogateModel         (          ); }
inline void                     TreedGaussianProcess::updateSurrogateModel      (void                      ) {        _root -> updateSurrogateModel      (          ); }
inline double                   TreedGaussianProcess::evaluateCriteria          (const vectord& query      ) { return _root -> evaluateCriteria          (query     ); }
inline void                     TreedGaussianProcess::updateCriteria            (const vectord& query      ) {        _root -> updateCriteria            (query     ); }
inline bool                     TreedGaussianProcess::criteriaRequiresComparison(void                      ) { return _root -> criteriaRequiresComparison(          ); }
inline void                     TreedGaussianProcess::setFirstCriterium         (void                      ) {        _root -> setFirstCriterium         (          ); }
inline bool                     TreedGaussianProcess::setNextCriterium          (const vectord& prevResult ) { return _root -> setNextCriterium          (prevResult); }
inline std::string              TreedGaussianProcess::getBestCriteria           (      vectord& best       ) { return _root -> getBestCriteria           (best      ); }
inline ProbabilityDistribution* TreedGaussianProcess::getPrediction             (const vectord& query      ) { return _root -> getPrediction             (query     ); }
inline uint                     TreedGaussianProcess::getNumberOfLeafs          (void                      ) { return _root -> getNumberOfLeafs          (          ); }
inline std::vector<TGPNode*>    TreedGaussianProcess::getLeafsPtrs              (void                      ) { std::vector<TGPNode*> result;
                                                                                                                      _root -> getLeafsPtrs              (result    );
                                                                                                               return result;                                          }

} // End: namespace bayesopt

#endif // _TREEDGAUSSIANPROCESS_H_
