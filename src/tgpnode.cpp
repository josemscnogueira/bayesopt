/**************************************************************************************************
 *  File:    tgpnode.cpp                                                                          *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/


/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include <cmath>

#include "log.hpp"
#include "ublas_trace.hpp"
#include "tgpnode.hpp"



/**************************************************************************************************
 *  Namespace: bayesopt                                                                           *
 **************************************************************************************************/
namespace bayesopt
{


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
TGPNode::TGPNode(TgpParameters& tgp_params, Parameters& params, randEngine& eng, Dataset& data, MeanModel& mean)
:
tgpparams (tgp_params   ),
parameters(params       ),
engine    (eng          )
{
    // Object Initialization
    _population    = new(uint);
    (*_population) = 1;
    _index         = 0;
    _type          = LEAF;

    // Bounds initialization
    _lower_bound = vectord(tgpparams.dimensions, 0.0);
    _upper_bound = vectord(tgpparams.dimensions, 1.0);

    // Update the node path to root
    _path.push_back(0);

    // Set parent to NULL
    _parent  = NULL;

    // Set childs to null
    _l_child = NULL;
    _r_child = NULL;

    // Create Dataset and MeanModel
    _data = new DatasetDist  (&data);
    _mean = new MeanModelDist(&mean);

    // Set wheights
    _wheights.clear    ( );
    _wheights.push_back(1);
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
TGPNode::TGPNode(TGPNode* parent, std::vector<uint>& owned, bool side, uint feature, double threshold)
:
tgpparams (parent -> tgpparams    ),
parameters(parent -> parameters   ),
engine    (parent -> engine       )
{
    // Initialize bounds
    _upper_bound = parent -> _upper_bound;
    _lower_bound = parent -> _lower_bound;

    _population  = parent -> _population;
    _index       = (*_population);
    _type        = LEAF;

    // Update bounds
    if (side == LEFT) _upper_bound(feature) = threshold;
    else              _lower_bound(feature) = threshold;

    // Update Path
    _path = parent -> _path;
    _path.push_back(_index);

    // Update Parent
    _parent  = parent;

    // Set childs to null
    _l_child = NULL;
    _r_child = NULL;

    // Create and update dataset and mean model
    _data = new DatasetDist  (parent -> _data, owned);
    _mean = new MeanModelDist(parent -> _mean, owned);

    // Increment population
    (*_population) += 1;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: Destructor                                                                       *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
TGPNode::~TGPNode(void)
{
    // Free Vectors
    _wheights    .clear();
    _outside_data.clear();
    _path        .clear();
    _lower_bound .clear();
    _upper_bound .clear();

    // Free bayesopt gps and and learning crit
    _gp  .clear();
    _crit.clear();

    // Free scoped_ptr
    _mcmc_sampler.reset();

    // Free tree in depth
    if (_l_child != NULL) delete _l_child; _l_child = NULL;
    if (_r_child != NULL) delete _r_child; _r_child = NULL;

    // Free data
    if (_mean != NULL ) delete _mean; _mean = NULL;
    if (_data != NULL ) delete _data; _data = NULL;

    if (_index == 0)
    {
        delete _population;
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getRoot                                                                          *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
TGPNode* TGPNode::getRoot(void)
{
    if (_parent == NULL) return this;
    else                 return _parent -> getRoot();
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getSplits                                                                        *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::getSplits(std::vector<double>& splits)
{
    if ( (_type == NON_LEAF) && (tgpparams.dimensions == 1) )
    {
        splits.push_back(_threshold);

        _l_child -> getSplits(splits);
        _r_child -> getSplits(splits);
    }
}

/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getNumberOfLeafs                                                                        *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
uint TGPNode::getNumberOfLeafs(void)
{
    uint result = 0;

    if (_type == LEAF)
    {
        return 1;
    }
    else // _type == NON_LEAF
    {
        result += _l_child -> getNumberOfLeafs();
        result += _r_child -> getNumberOfLeafs();
    }

    return result;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getLeafsPtrs                                                                     *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::getLeafsPtrs(std::vector<TGPNode*>& result)
{
    if (_type == LEAF)
    {
        result.push_back(this);
    }
    else // _type == NON_LEAF
    {
        _l_child -> getLeafsPtrs(result);
        _r_child -> getLeafsPtrs(result);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getBounds                                                                        *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::getBounds(vectord& lower, vectord& upper)
{
    lower = _lower_bound;
    upper = _upper_bound;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getHyperParameters                                                               *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::getHyperParameters(std::vector<vectord>& lower_bounds, std::vector<vectord>& upper_bounds, std::vector< std::vector<vectord> >& hypers)
{
    if (_index == 0)
    {
        lower_bounds.clear();
        upper_bounds.clear();
        hypers      .clear();
    }

    if (_type == LEAF)
    {
        std::vector<vectord> hyper_particles;
                             hyper_particles.clear();

        // Fill particle vector
        for (uint index = 0; index < tgpparams.mcmc_particles; index += 1)
        {
            hyper_particles.push_back(_gp[index].getHyperParameters());
        }

        // Update outputs
        hypers      .push_back(hyper_particles);
        lower_bounds.push_back(_lower_bound   );
        upper_bounds.push_back(_upper_bound   );
    }
    else // _type == NON_LEAF
    {
        _l_child -> getHyperParameters(lower_bounds, upper_bounds, hypers);
        _r_child -> getHyperParameters(lower_bounds, upper_bounds, hypers);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: loadHyperParameters                                                              *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::loadHyperParameters(std::vector<vectord>& lower_bounds, std::vector<vectord>& upper_bounds, std::vector< std::vector<vectord> >& hypers)
{
    if (_type == LEAF)
    {
        uint                hyperparameter_size = _gp[0].nHyperParameters();
        std::vector<double> superposition       = calculateSuperposition(lower_bounds, upper_bounds);

        for (uint particle = 0; particle < tgpparams.mcmc_particles; particle += 1)
        {
            vectord hyperparameters = vectord(hyperparameter_size, 0.0);

            for (uint leaf = 0; leaf < hypers.size(); leaf += 1)
            {
                hyperparameters += superposition[leaf] * (hypers[leaf])[particle];
            }

            _gp[particle].setHyperParameters(hyperparameters);
        }
    }
    else
    {
        _l_child -> loadHyperParameters(lower_bounds, upper_bounds, hypers);
        _r_child -> loadHyperParameters(lower_bounds, upper_bounds, hypers);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: copyHyperparameters                                                              *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::copyHyperparameters(TGPNode* old_root, TGPNode* new_root)
{
    std::vector< vectord              > lower_bounds;
    std::vector< vectord              > upper_bounds;
    std::vector< std::vector<vectord> > hypers;

    old_root -> getHyperParameters (lower_bounds, upper_bounds, hypers);
    new_root -> loadHyperParameters(lower_bounds, upper_bounds, hypers);
    new_root -> fitSurrogateModel  (                                  );
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: calculateSuperposition                                                           *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
std::vector<double> TGPNode::calculateSuperposition(std::vector<vectord>& lower_bounds, std::vector<vectord>& upper_bounds)
{
    std::vector<double> result;
                        result.clear();

    for (uint leaf = 0; leaf < lower_bounds.size(); leaf += 1)
    {
        double superposition = 1;

        for (uint dimension = 0; dimension < _lower_bound.size(); dimension += 1)
        {
            double upper_min = std::min(upper_bounds[leaf](dimension), _upper_bound(dimension));
            double lower_max = std::max(lower_bounds[leaf](dimension), _lower_bound(dimension));

            if (lower_max >= upper_min) { superposition = 0.0; break; }

            superposition *= ( (upper_min - lower_max) / (_upper_bound(dimension) - _lower_bound(dimension)) );
        }

        result.push_back(superposition);
    }

    return result;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setTree                                                                          *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::setTree(void)
{
    if (_type == LEAF)
    {
        // Update max and min values
        _data -> setMinMax();

        // Set Surrogate models and Criteria
        setSurrogateModel();
        setCriteria();

        // Calculate Wheghts and set pointers to leaf data
        calculateWheights(1);

        // Set MCMCSampler
        _mcmc_sampler   .reset        (new MCMCSampler(this, _gp[0].nHyperParameters(), engine));
        _mcmc_sampler -> setNParticles(tgpparams.mcmc_particles);
        _mcmc_sampler -> setNBurnOut  (100);
    }
    else
    {
        _l_child -> setTree();
        _r_child -> setTree();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: isEqual                                                                          *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
bool TGPNode::isEqual(TGPNode* node)
{
    if (_index == 0)
    {
        if ((*_population) != (*(node -> _population))) return false;
    }

    if (_type == NON_LEAF)
    {
        // Verify if split is equal
        if ( ( _feature   != node -> _feature   ) ||
             ( _threshold != node -> _threshold )   )
        {
            return false;
        }

        // Verify if left  child exists or not for both
        if ( ( ( _l_child == NULL ) && ( node -> _l_child != NULL ) ) ||
             ( ( _l_child != NULL ) && ( node -> _l_child == NULL ) )    )
        {
            return false;
        }


        // Verify if right child exists or not for both
        if ( ( ( _r_child == NULL ) && ( node -> _r_child != NULL ) ) ||
             ( ( _r_child != NULL ) && ( node -> _r_child == NULL ) )    )
        {
            return false;
        }

        // Check if left  child is equal
        if ( ( _l_child != NULL ) && ( _r_child != NULL ) )
        {
            if ( _l_child -> isEqual(node -> _l_child) == false)
            {
                return false;
            }
        }

        // Check if right child is equal
        if ( ( _r_child != NULL ) && ( _r_child != NULL ) )
        {
            if ( _r_child -> isEqual(node -> _r_child) == false)
            {
                return false;
            }
        }
    }

    // Return true if everything checks out (means this is a leaf)
    return true;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: addSample                                                                        *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::addSample(const vectord &x, double y)
{
    // If the node is a leaf, add the point
    if (_type == LEAF)
    {
        if (isOutOfBounds(x))
        {
            FILE_LOG(logERROR) << "Error in TGPNode::addSample: Sample is out of bounds.";
            exit(-1);
        }

        // Add data
        _data -> addSample  (x,y);
        _mean -> addNewPoint(x  );
    }
    // If not, pass the point down the tree
    else
    {
        if (x(_feature) < _threshold)  _l_child -> addSample(x,y);
        else                           _r_child -> addSample(x,y);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: isOutOfBounds                                                                    *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
bool TGPNode::isOutOfBounds(const vectord& x)
{
    // Check if query satisfies bounds of the node
    for (uint index = 0; index < tgpparams.dimensions; index += 1)
    {
        if ( x(index) > _upper_bound[index] ||
             x(index) < _lower_bound[index]   )
        {
            return true;
        }
    }

    return false;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: createSplit                                                                      *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::createSplit(std::vector<uint>& div_left, std::vector<uint>& div_right, uint feature, double threshold)
{
    // Create child nodes
    if (_l_child != NULL) delete _l_child; _l_child = NULL;
    if (_r_child != NULL) delete _r_child; _r_child = NULL;

    _l_child = new TGPNode(this, div_left , LEFT , feature, threshold);
    _r_child = new TGPNode(this, div_right, RIGHT, feature, threshold);

    // Update to non leaf node
    _type      = NON_LEAF;
    _threshold = threshold;
    _feature   = feature;

    // Delete Node data
    delete _data; _data = NULL;
    delete _mean; _mean = NULL;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: calculateUncertaintyReduction                                                     *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
double TGPNode::calculateUncertaintyReduction(std::vector<uint>& div_left, std::vector<uint>& div_right, uint feature, double threshold)
{
    // // Calculate the uncertainty of the current node
    double result, gamma0, gamma1, gamma2;

    // Calculate the uncertainty of the current node
    gamma0 = _data -> calculateUncertainty();

    // // Calculate the reduction of the split
    gamma1 = _data -> calculateUncertainty(div_left );
    gamma2 = _data -> calculateUncertainty(div_right);

    result = gamma0 - gamma1 - gamma2;

    // Return value
    return result;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: reduceNodeUncertainty                                                             *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
double TGPNode::reduceNodeUncertainty(void)
{
    // If not leaf, no reduction
    if (_type != LEAF) return 0.0;

    // If there is no sufficient data, no reduction
    if (_data -> getNSamples() < (tgpparams.min_data_per_leaf * 2)) return 0.0;

    // Search Variables
    double            uncertainty_reduction = 0.0;
    bool              invalid_split;
    std::vector<uint> division_left;
    std::vector<uint> division_right;
    std::vector<uint> best_division_left;
    std::vector<uint> best_division_right;
    uint              best_feature;
    double            best_threshold;

    for (uint sample = 0; sample < _data -> getNSamples(); sample  += 1)
    {
        vectord split = _data -> getSampleX(sample);

        for (uint feature = 0; feature < tgpparams.dimensions; feature += 1)
        {
            // Don't consider leafs with small size
            double input_noise = std::sqrt(parameters.crit_params[4 + (feature * tgpparams.dimensions) + feature]);

            if (input_noise == 0) input_noise = 0.01;

            if ( ((       split(feature) - _lower_bound[feature]) < input_noise) ||
                 ((_upper_bound[feature] -        split(feature)) < input_noise)   )
            {
                continue;
            }

            // If it is big enough, start the split
            division_left .clear();
            division_right.clear();
            invalid_split  = false;

            for (uint index = 0; index < _data -> getNSamples(); index += 1)
            {
                // Split according to a feature and a specific point
                if (index != sample)
                {
                    if      ( _data -> getSampleX(index)(feature) < split(feature))
                    {
                        // Left split
                        division_left .push_back(_data -> getOwnedData(index));
                    }
                    else if ( _data -> getSampleX(index)(feature) > split(feature))
                    {
                        // Right split
                        division_right.push_back(_data -> getOwnedData(index));
                    }
                    else
                    {
                        // There's an invalid splot for this feature
                        invalid_split = true;
                        break;
                    }
                }
                else
                {
                    division_left .push_back(_data -> getOwnedData(index));
                    division_right.push_back(_data -> getOwnedData(index));
                }
            }

            // If the split is invalid, proceed to next feature
            if (invalid_split) continue;

            // Not enough data to split
            if ( division_left .size() < tgpparams.min_data_per_leaf ||
                 division_right.size() < tgpparams.min_data_per_leaf   ) continue;

            // If there is enought data to split
            double uncertainty = calculateUncertaintyReduction(division_left, division_right, feature, split(feature));

            // The split isn't doesn't minimezes uncertainty reduction
            if (uncertainty > uncertainty_reduction)
            {
                uncertainty_reduction = uncertainty;

                best_division_left  = division_left;
                best_division_right = division_right;
                best_feature        = feature;
                best_threshold      = split(feature);
            }
        }
    }

    // Evaluate if there was an uncertainty reduction
    if (uncertainty_reduction > 0.0)
    {
        // Split Node
        createSplit(best_division_left, best_division_right, best_feature, best_threshold);

        return uncertainty_reduction;
    }
    else
    {
        return 0.0;
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: reduceTreeUncertainty                                                            *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::reduceTreeUncertainty(void)
{
    if(_type == LEAF)
    {
        if (reduceNodeUncertainty() > 0.0)
        {
                _l_child -> reduceTreeUncertainty();
                _r_child -> reduceTreeUncertainty();
        }
    }
    else
    {
        _l_child -> reduceTreeUncertainty();
        _r_child -> reduceTreeUncertainty();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setSurrogateModel                                                                *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::setSurrogateModel(void)
{
    if (_type == LEAF)
    {
        _gp.clear();

        for (uint index = 0; index < tgpparams.mcmc_particles; index += 1)
        {
            _gp.push_back(new GaussianProcess(tgpparams.dimensions, parameters, *_data, *_mean, engine));
        }

    }
    else // _type == NON_LEAF
    {
        _l_child -> setSurrogateModel();
        _r_child -> setSurrogateModel();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: fitSurrogateModel                                                                *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::fitSurrogateModel(void)
{
    if (_type == LEAF)
    {
        if (_data -> getNSamples() == 0) return;

        for (uint index = 0; index <  _gp.size(); index += 1)
        {
            _gp[index].fitSurrogateModel();
        }
    }
    else
    {
        _l_child -> fitSurrogateModel();
        _r_child -> fitSurrogateModel();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: updateSurrogateModel                                                             *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::updateSurrogateModel(void)
{
    if (_type == LEAF)
    {
        for (GPVect::iterator gp = _gp.begin(); gp != _gp.end(); gp += 1)
        {
            gp -> updateSurrogateModel();
        }
    }
    else
    {
        _l_child -> updateSurrogateModel();
        _r_child -> updateSurrogateModel();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setCriteria                                                                      *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::setCriteria(void)
{
    if (_type == LEAF)
    {
        CriteriaFactory mcfactory;

        _crit.clear();

        for (uint index = 0; index < tgpparams.mcmc_particles; index += 1)
        {
            _crit.push_back(mcfactory.create(parameters.crit_name, &_gp[index]));
            _crit[index].setRandomEngine(engine);

            if (_crit[index].nParameters() == parameters.crit_params.size())
            {
                _crit[index].setParameters(parameters.crit_params);
            }
            else
            {
                if (parameters.crit_params.size() != 0)
                {
                    FILE_LOG(logERROR) << "Expected " << _crit[index].nParameters()
                                       << " parameters. Got "
                                       << parameters.crit_params.size() << " instead.";
                }

                FILE_LOG(logERROR)     << "Using default parameters for criteria.";

                exit(-1);
            }
        }
    }
    else // _type == NON_LEAF
    {
        _l_child -> setCriteria();
        _r_child -> setCriteria();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: evaluateCriteria                                                                 *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
double TGPNode::evaluateCriteria(const vectord& query)
{
    if (_type == LEAF)
    {
        if (isOutOfBounds(query))
        {
            FILE_LOG(logERROR) << "Error in TGPNode::evaluateCriteria: Query is out of bounds.";
            exit(-1);
        }

        double result = 0.0;

        for (uint index = 0; index < _crit.size(); index += 1)
        {
            result += _crit[index].evaluate(query);
        }

        return ( result / static_cast<double>(tgpparams.mcmc_particles) );
    }
    else // _type == NON_LEAF
    {
        if ( query(_feature) < _threshold ) return _l_child -> evaluateCriteria(query);
        else                                return _r_child -> evaluateCriteria(query);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: updateCriteria                                                                   *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::updateCriteria(const vectord& query)
{
    if (_type == LEAF)
    {
        if (isOutOfBounds(query))
        {
            FILE_LOG(logERROR) << "Error in TGPNode::updateCriteria: Query is out of bounds.";
            exit(-1);
        }

        for (CritVect::iterator crit = _crit.begin(); crit != _crit.end(); crit += 1)
        {
            crit -> update(query);
        }
    }
    else // _type == NON_LEAF
    {
        if ( query(_feature) < _threshold ) _l_child -> updateCriteria(query);
        else                                _r_child -> updateCriteria(query);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: criteriaRequiresComparison                                                       *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
bool TGPNode::criteriaRequiresComparison(void)
{
    if (_type == LEAF) return _crit[0].requireComparison();
    else               return _l_child -> criteriaRequiresComparison();
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setFirstCriterium                                                                *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::setFirstCriterium(void)
{
    if (_type == LEAF)
    {
        for (CritVect::iterator crit = _crit.begin(); crit != _crit.end(); crit += 1)
        {
            crit -> initialCriteria();
        }
    }
    else // _type == NON_LEAF
    {
        _l_child -> setFirstCriterium();
        _r_child -> setFirstCriterium();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setNextCriterium                                                                 *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
bool TGPNode::setNextCriterium(const vectord& prevResult)
{
    bool rotated;

    pushResult    (prevResult);
    rotateCriteria(rotated   );

    return rotated;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: pushResult                                                                       *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::pushResult(const vectord& prevResult)
{
    if (_type == LEAF)
    {
        if (isOutOfBounds(prevResult))
        {
            FILE_LOG(logERROR) << "Error in TGPNode::pushResult: Query is out of bounds.";
            exit(-1);
        }

        _crit[0].pushResult(prevResult);
    }
    else // _type == NON_LEAF
    {
        if ( prevResult(_feature) < _threshold ) _l_child -> pushResult(prevResult);
        else                                     _r_child -> pushResult(prevResult);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: rotateCriteria                                                                   *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::rotateCriteria(bool& rotated)
{
    if (_type == LEAF)
    {
        for (CritVect::iterator crit = _crit.begin(); crit != _crit.end(); crit += 1)
        {
            rotated = crit -> rotateCriteria();
        }
    }
    else
    {
        _l_child -> rotateCriteria(rotated);
        _r_child -> rotateCriteria(rotated);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getBestCriteria                                                                  *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
std::string TGPNode::getBestCriteria(vectord& best)
{
    if (_type == LEAF)
    {
        if (isOutOfBounds(best))
        {
            FILE_LOG(logERROR) << "Error in TGPNode::getBestCriteria: Query is out of bounds.";
            exit(-1);
        }

        return _crit[0].getBestCriteria(best);
    }
    else // _type == NON_LEAF
    {
        if ( best(_feature) < _threshold ) return _l_child -> getBestCriteria(best);
        else                               return _r_child -> getBestCriteria(best);
    }

    return "Never reached";
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getPrediction                                                                    *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
ProbabilityDistribution* TGPNode::getPrediction(const vectord& query)
{
    if (_type == LEAF)
    {
        return _gp[0].prediction(query);
    }
    else // _type == NON_LEAF
    {
        if ( query(_feature) < _threshold ) return _l_child -> getPrediction(query);
        else                                return _r_child -> getPrediction(query);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: setHyperParameters                                                               *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::setHyperParameters(const vectord &theta)
{
    if (_type == LEAF)
    {
        _gp[0].setHyperParameters(theta);
    }
    else // _type == NON_LEAF
    {
        exit(-1);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: calculateWheights                                                                *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::calculateWheights(uint mode)
{
    if (mode == 0)
    {
        TGPNode* root = getRoot();

        root -> calculateWheights(1);
    }
    else
    {
        if (_type == LEAF)
        {
            // Clear Previous Wheights
            _wheights    .clear();
            _outside_data.clear();

            // Set first node as current first node in the path
            TGPNode* node        = this;
            uint     child_index = UINT_MAX;

            // Calculate wheights along the path
            while (node != NULL)
            {
                // Calculate wheight
                _wheights.push_back(calculateWheight(node));

                // Update Other nodes' data
                _outside_data.push_back(std::vector<std::vector<uint>*>());

                node -> getOutsideData(child_index, _outside_data.back());

                // Next node
                child_index = node -> _index;
                node        = node -> _parent;
            }
        }
        else // _type == NON_LEAF
        {
            _l_child -> calculateWheights(1);
            _r_child -> calculateWheights(1);
        }
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: calculateWheight                                                                 *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
double TGPNode::calculateWheight(TGPNode* node)
{
    double result    = _path.size(); // result = depth_i + 1
    uint   index     = 0;

    while ( (index <         _path.size()) &&
            (index < node -> _path.size())    )
    {
        if ( _path[index] == node -> _path[index] ) result -= 1; // depth_ij += 1

        index += 1;
    }

    result = 2.0 / (result + 1); // result = depth_i + 1 - (depth_j + 1), den = result + 1

    result = std::pow(result, tgpparams.wheight_power);

    return result;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: getOutsideData                                                                   *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::getOutsideData(uint index, std::vector<std::vector<uint>*>& data)
{
    // If node was already visited for other wheight, just return
    if (_index == index) return;

    // Else, node was not yet included
    if (_type == LEAF)
    {
        data.push_back(_data -> getOwnedDataPtr());
    }
    else // _type == NON_LEAF
    {
        _l_child -> getOutsideData(index, data);
        _r_child -> getOutsideData(index, data);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: evaluate                                                                         *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
double TGPNode::evaluate(const vectord& x)
{
    if (_type == LEAF)
    {
        setHyperParameters(x);

        return evaluateKernelParams();
    }
    else // _type == NON_LEAF
    {
        exit(-1);
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: evaluateKernelParams                                                             *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
double TGPNode::evaluateKernelParams(void)
{
    // Calculate kernel prior
    double  result  = 0.0;
            result += - _gp[0].getKernelPrior(); // TODO: !!!!

	// Add Log Marginal Likelihoods
	if (tgpparams.wheight_power > 3)
	{
		result += negativeLogLikelihood(0);
	}
	else
	{
		for (uint index = 0; index < _wheights.size(); index += 1)
		{
			result += _wheights[index] * negativeLogLikelihood(index);
		}
	}

    return result;
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: negativeLogLikelihood                                                            *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
double TGPNode::negativeLogLikelihood(uint index)
{
    namespace ublas = boost::numeric::ublas;
    // According with algorithm 2.1 from Gaussian Processes for Machine Learning, Rasmussen

    DatasetDist*   data = new DatasetDist  (_data, _outside_data[index]     );
    MeanModelDist* mean = new MeanModelDist(_mean, data -> getOwnedDataPtr());

    if ( (data == NULL) || (mean == NULL) ) exit(-1);

    const double  sigma = _gp[0].getSignalVariance();
    const matrixd K     = computeCorrMatrix(*data);
    const size_t  n     = K.size1();
          matrixd L(n,n);

    // FILE_LOG(logINFO) << "K" << std::endl << K << std::endl;

    // Compute Cholesky Matrix L
    utils::cholesky_decompose(K,L);

    // Calculate alpha
    vectord alpha(data -> getSamplesY() - mean -> muTimesFeat());

    inplace_solve(L, alpha, ublas::lower_tag());

    // Return loglikelihood
    double result  = ublas::inner_prod(alpha,alpha)/(2*sigma); // TODO: VERIFY this
           result += utils::log_trace(L);
           result += n * std::log(2*M_PI) / 2;

    delete data;
    delete mean;

    return result;
}

double TGPNode::mgpdf(const vectord& x, const vectord& mu, const matrixd& S, double p)
{
    size_t  n = mu.size();
    matrixd L(n,n);

    bayesopt::utils::cholesky_decompose(S,L);

    vectord u = mu - x;

    inplace_solve(L,u,boost::numeric::ublas::lower_tag());

    double quad = boost::numeric::ublas::inner_prod(u,u) /2.0;
    double det  = bayesopt::utils::log_trace(L);
    double cst  = n * std::log(M_PI) / 2.0;

    return p * std::exp(-quad-det-cst);
}



/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: updateHyperParameters                                                            *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TGPNode::updateHyperParameters(void)
{
    if (_type == LEAF)
    {
        if (_data -> getNSamples() == 0) return;

        size_t  last      = _gp.size() - 1;
        vectord lastTheta = _gp[last].getHyperParameters();

        FILE_LOG(logDEBUG) << "Hyperparameter Optimization -> Leaf with index = " << _index << "; and visible nodes = " << _wheights.size();

        _mcmc_sampler -> run(lastTheta);

        for(size_t index = 0; index < tgpparams.mcmc_particles; index += 1)
        {
            _gp[index].setHyperParameters(_mcmc_sampler -> getParticle(index));
        }

    }
    else // _type == NON_LEAF
    {
        _l_child -> updateHyperParameters();
        _r_child -> updateHyperParameters();
    }
}


/**************************************************************************************************
 *  Procedure                                                                                     *
 *                                                                                                *
 *  Description: print                                                                            *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
std::string whiteSpacesString(uint tab_counter)
{
    std::string s;

    if (tab_counter == 0) return s;

    for(uint index = 0; index < tab_counter; index += 1)
    {
        s += "    ";
    }

    return s;
}


void TGPNode::print(void)
{
    // Node Attributes
    FILE_LOG(logINFO) << whiteSpacesString(_path.size() - 1) << "{";
    FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Index     = " <<  _index            << ",";
    FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Type      = " <<  _type             << ",";
    FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Depth     = " << (_path.size() - 1) << ",";

    if( _type == LEAF)
    {
        // Data Attributes
        if (_data)
        {
            FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Data      = " << _data -> getNSamples();
            FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) << "{"                                          ;

            for(uint index = 0; index < _data -> getNSamples(); index += 1)
            {
                FILE_LOG(logINFO) << whiteSpacesString(_path.size() + 1) << _data -> getSampleX(index) << ", " << _data -> getSampleY(index);
            }

            FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) << "}";

            FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Wheights  = " << _wheights.size();
            FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) << "{"                                    ;

            for(uint index = 0; index < _wheights.size(); index += 1)
            {
                FILE_LOG(logINFO) << whiteSpacesString(_path.size() + 1) << _wheights[index] << ", ";
            }
            FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) << "}";
                                  ;
            for (uint index1 = 0; index1 < _outside_data.size(); index1 += 1)
            {
                FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Outside Data  = " << index1;
                FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) << "{";

                DatasetDist*   data = new DatasetDist  (_data, _outside_data[index1]);
                vecOfvec xx =  data -> getSamplesX();

                for (uint sample = 0; sample < xx.size(); sample += 1)
                {
                    FILE_LOG(logINFO) << whiteSpacesString(_path.size() + 1 ) << xx[sample] << ", ";
                }

                delete data;

                FILE_LOG(logINFO) << whiteSpacesString(_path.size() + 0 ) << "}";
            }
        }
    }
    else
    {
        FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Feature   = " <<  _feature          << ",";
        FILE_LOG(logINFO) << whiteSpacesString(_path.size()    ) <<     " Threshold = " <<  _threshold        << ",";

        // Child print calls
        if (_l_child != NULL) _l_child -> print();
        if (_r_child != NULL) _r_child -> print();
    }

    // End
    FILE_LOG(logINFO) << whiteSpacesString(_path.size() - 1) << "}";
}



} // END of namespace bayesopt
