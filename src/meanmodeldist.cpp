/**************************************************************************************************
 *  File:    meanmodeldist.cpp                                                                    *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/


/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include "meanmodeldist.hpp"
#include "log.hpp"


/**************************************************************************************************
 *  Namespace: bayesopt                                                                           *
 **************************************************************************************************/
namespace bayesopt
{


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : MeanModelDist                                                                    *
 **************************************************************************************************/
MeanModelDist::MeanModelDist(MeanModel* mean)
{
    // Update all data
    _alldata = mean;

    // Own all data
    if (mean -> getNSamples() > 0)
    {
        for (uint index = 0; index < mean -> getNSamples(); index += 1)
        {
            _owned.push_back(index);
        }
    }
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : MeanModelDist                                                                    *
 **************************************************************************************************/
MeanModelDist::MeanModelDist(MeanModelDist* mean, std::vector<uint>& owned)
{
    // Update all data
    _alldata = mean -> _alldata;

    // Update owned data
    _owned = owned;
}

/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : MeanModelDist                                                                    *
 **************************************************************************************************/
MeanModelDist::MeanModelDist(MeanModelDist* mean, std::vector<uint>* owned)
{
    // Update all data
    _alldata = mean -> _alldata;

    // Update owned data
    _owned = (*owned);
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: setPoints                                                                        *
 *  Class      : MeanModelDist                                                                    *
 **************************************************************************************************/
void MeanModelDist::setPoints(const vecOfvec& x)
{
    // Update all data
    _alldata -> setPoints(x);

    // Update owned data
    _owned.clear();

    for(uint index = 0; index < x.size(); index += 1)
    {
        _owned.push_back(index);
    }
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: addNewPoint                                                                      *
 *  Class      : MeanModelDist                                                                    *
 **************************************************************************************************/
void MeanModelDist::addNewPoint(const vectord &x)
{
    // Update all data
    _alldata -> addNewPoint(x);

    // Update owned data
    _owned.push_back( _alldata -> mFeatM.size2() - 1 );
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: getmFeatM                                                                        *
 *  Class      : MeanModelDist                                                                    *
 **************************************************************************************************/
matrixd MeanModelDist::getmFeatM(void)
{
    using boost::numeric::ublas::column;

    matrixd result(_alldata -> mFeatM.size1(), _owned.size());

    for (uint index = 0; index < _owned.size(); index += 1)
    {
        column(result, index) = column(_alldata -> mFeatM, _owned[index]);
    }

    return result;
}

} // END of namespace bayesopt
