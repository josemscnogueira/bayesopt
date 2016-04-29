/**************************************************************************************************
 *  File:    tgpnode.h                                                                            *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/

#ifndef __TGPNODE_HPP__
#define __TGPNODE_HPP__

#include "ublas_cholesky.hpp"
#include "ublas_trace.hpp"

/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include <vector>                              // Global
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>

#include "mcmc_sampler.hpp"                    // Bayesopt
#include "gaussian_process.hpp"
#include "criteria_functors.hpp"
#include "datasetdist.hpp"
#include "meanmodeldist.hpp"
#include "tgpparameters.h"


/**************************************************************************************************
 *  Namespace: bayesopt                                                                           *
 **************************************************************************************************/
namespace bayesopt
{


/**************************************************************************************************
 *  Enumerations                                                                                  *
 **************************************************************************************************/
enum
{
    LEAF     = 0,
    NON_LEAF = 1,
};

enum
{
    LEFT     = false,
    RIGHT    = true ,
};


/**************************************************************************************************
 *  Class: TGPNode                                                                                *
 **************************************************************************************************/
class TGPNode : public RBOptimizable
{
public:
    // Typedefs
    typedef boost::scoped_ptr<MCMCSampler>           MCMCSamplerPtr;
    typedef boost::ptr_vector<GaussianProcess>       GPVect;
    typedef boost::ptr_vector<Criteria>              CritVect;


    // Constructors
    TGPNode(tgp_parameters& tgp_params, bopt_params& params      , randEngine& eng , Dataset& data   , MeanModel& mean);
    TGPNode(TGPNode* parent           , std::vector<uint>& owned , bool        side, uint     feature, double     threshold);

    // Destructor
    ~TGPNode(void);

    // Static Methods
    static void              copyHyperparameters       (TGPNode* old_root, TGPNode* new_root);

    //Methods
    TGPNode*                 getRoot                   (void);
    void                     getSplits                 (std::vector<double>&   splits);
    uint                     getNumberOfLeafs          (void);
    void                     getLeafsPtrs              (std::vector<TGPNode*>& result);
    void                     setTree                   (void);
    bool                     isEqual                   (TGPNode* node);
    void                     reduceTreeUncertainty     (void);
    void                     getBounds                 (vectord& lower, vectord& upper);

    void                     addSample                 (const vectord &x, double y);

    void                     fitSurrogateModel         (void);
    void                     updateSurrogateModel      (void);
    void                     updateHyperParameters     (void);

    void                     getHyperParameters        (std::vector<vectord>& lower_bounds, std::vector<vectord>& upper_bounds, std::vector< std::vector<vectord> >& hypers);
    void                     loadHyperParameters       (std::vector<vectord>& lower_bounds, std::vector<vectord>& upper_bounds, std::vector< std::vector<vectord> >& hypers);

    void                     setCriteria               (void);
    double                   evaluateCriteria          (const vectord& query);
    void                     updateCriteria            (const vectord& query);
    bool                     criteriaRequiresComparison(void);
    void                     setFirstCriterium         (void);
    bool                     setNextCriterium          (const vectord& prevResult);
    std::string              getBestCriteria           (      vectord& best);

    ProbabilityDistribution* getPrediction             (const vectord& query);

    double                   evaluate                  (const vectord& x    );

    void                     print                     (void);

    static double            mgpdf(const vectord& x, const vectord& mu, const matrixd& S, double p);


protected:
    // Attributes
    DatasetDist*                                      _data;
    MeanModelDist*                                    _mean;

    std::vector<double>                               _wheights;
    std::vector< std::vector< std::vector<uint>* > >  _outside_data;

private:
    // Object Referenfes
    tgp_parameters&                                   tgpparams;
    bopt_params&                                      parameters;
    randEngine&                                       engine;

    // Node Attributes
    uint*                                             _population;
    uint                                              _index;
    uint                                              _type;
    uint                                              _feature;
    double                                            _threshold;
    std::vector<uint>                                 _path;

    vectord                                           _lower_bound;
    vectord                                           _upper_bound;

    GPVect                                            _gp;
    CritVect                                          _crit;

    TGPNode*                                          _parent;
    TGPNode*                                          _l_child;
    TGPNode*                                          _r_child;

    MCMCSamplerPtr                                    _mcmc_sampler;

    // Methods
    bool                isOutOfBounds                (const vectord& x);
    void                createSplit                  (std::vector<uint>&  div_left, std::vector<uint>& div_right, uint feature, double threshold);
    double              calculateUncertainty         (std::vector<uint>&  data);
    double              calculateUncertainty         (void);
    double              calculateUncertaintyReduction(std::vector<uint>&  div_left, std::vector<uint>& div_right, uint feature, double threshold);
    double              reduceNodeUncertainty        (void);
    void                setSurrogateModel            (void);

    // Criteria methods
    void                pushResult                   (const vectord& prevResult);
    void                rotateCriteria               (      bool&    rotated);

    // Hyperparameter optimization methods
    void                setHyperParameters           (const vectord& theta);
    double              evaluateKernelParams         (void                );
    double              negativeLogLikelihood        (uint           index);
    matrixd             computeCorrMatrix            (DatasetDist&   data );

    // Leaf Wheights
    void                calculateWheights            (uint           mode = 0);
    double              calculateWheight             (TGPNode*       node    );
    void                getOutsideData               (uint index             , std::vector<std::vector<uint>*>& data);

    // Other
    std::vector<double> calculateSuperposition       (std::vector<vectord>& lower_bounds, std::vector<vectord>& upper_bounds);

};

/**************************************************************************************************
 *  Inline methods : TGPNode                                                                      *
 **************************************************************************************************/
inline double  TGPNode::calculateUncertainty(std::vector<uint>& data) { return _data -> calculateUncertainty(data              ); }
inline double  TGPNode::calculateUncertainty(void              )      { return _data -> calculateUncertainty(                  ); }
inline matrixd TGPNode::computeCorrMatrix   (DatasetDist&  data)      { return _gp[0]  .calculateCorrMatrix(data.getSamplesX() ); }

} // End: namespace bayesopt

#endif // _TGPNode_h_
