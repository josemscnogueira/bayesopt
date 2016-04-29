/**************************************************************************************************
 *  File:    treedgaussianprocess.cpp                                                             *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/


/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include "treedgaussianprocess.hpp"


/**************************************************************************************************
 *  Namespace: bayesopt                                                                           *
 **************************************************************************************************/
namespace bayesopt
{

/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
TreedGaussianProcess::TreedGaussianProcess(TgpParameters& tgp_params, Parameters& params, randEngine& eng)
:
PosteriorModel(tgp_params.dimensions,params,eng),
tgpparams     (tgp_params                      ),
engine        (eng                             )
{
    // Create root
    _root = new TGPNode(tgpparams, params, eng, mData, mMean);
    _root -> setTree();
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Destructor                                                                       *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
TreedGaussianProcess::~TreedGaussianProcess()
{
    delete _root;
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: setSamples                                                                       *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TreedGaussianProcess::setSamples(const matrixd &x, const vectord &y)
{
    for(uint index = 0; index < x.size1(); index += 1)
    {
        _root -> addSample( row(x,index), y(index));
    }
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: setSample                                                                        *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TreedGaussianProcess::setSample(const vectord &x, double y)
{
    _root -> addSample(x,y);
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: setSample                                                                        *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
std::vector<double> TreedGaussianProcess::getSplits(void)
{
    std::vector<double> result;

    if (mDims != 1) return result;

    _root -> getSplits(result);

    return result;
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: updateHyperParameters                                                            *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TreedGaussianProcess::updateHyperParameters(void)
{
    TGPNode* new_root = new TGPNode(tgpparams, mParameters, engine, mData, mMean);
             new_root -> reduceTreeUncertainty();

    if ( _root -> isEqual(new_root) )
    {
        delete new_root;
    }
    else
    {
        // Create new TGP
        new_root -> setTree();

        // Load previous tree's hyperparameters
        TGPNode::copyHyperparameters(_root, new_root);

        // Free memory
        delete _root;

        // New tree assign
        _root = new_root;
    }

    _root -> updateHyperParameters();
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: print                                                                            *
 *  Class      : TGPNode                                                                          *
 **************************************************************************************************/
void TreedGaussianProcess::print(void)
{
    FILE_LOG(logINFO) << "/**************************************************************************************************" ;
    FILE_LOG(logINFO) << " *  TreedGaussianProcess                                                                          *" ;
    FILE_LOG(logINFO) << " **************************************************************************************************/";

    _root -> print();

    if (mData.getNSamples() > 0)
    {
        FILE_LOG(logINFO) << "/**************************************************************************************************";
        FILE_LOG(logINFO) << "Minimizer = " << getPointAtMinimum();
        FILE_LOG(logINFO) << "Minimum   = " << getValueAtMinimum();
    }

    FILE_LOG(logINFO) << "/**************************************************************************************************";
    FILE_LOG(logINFO) << " *                                                                                          [end] *";
    FILE_LOG(logINFO) << " **************************************************************************************************/";
}


} // END of namespace bayesopt
