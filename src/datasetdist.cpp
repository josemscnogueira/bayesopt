/**************************************************************************************************
 *  File:    datasetdist.cpp                                                                      *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/


/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include "datasetdist.hpp"



/**************************************************************************************************
 *  Namespace: bayesopt                                                                           *
 **************************************************************************************************/
namespace bayesopt
{


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
DatasetDist::DatasetDist(Dataset* data)
{
    // Update all data
    _alldata = data;

    // Own all data
    if (data -> getNSamples() > 0)
    {
        for (uint index = 0; index < data -> getNSamples(); index += 1)
        {
            _owned.push_back(index);
        }
    }

    // Update maximum and minimum indexes
    _min = 0;
    _max = 0;
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
DatasetDist::DatasetDist(DatasetDist* data, std::vector<uint>& owned)
{
    // Update all data
    _alldata = data -> _alldata;

    // Update owned vector
    _owned = owned;

    // Update Min and Max indexes
    _min = 0;
    _max = 0;
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: Constructor                                                                      *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
DatasetDist::DatasetDist(DatasetDist* data, std::vector<std::vector<uint>*> owned)
{
    // Update all data
    _alldata = data -> _alldata;

    // Update owned vector
    _owned.clear();

    // Reserve owned
    uint size = 0;

    for (uint index = 0; index < owned.size(); index += 1)
    {
        size += owned[index] -> size();
    }

    _owned.reserve(size);

    // Fill owned
    for (uint index = 0; index < owned.size(); index += 1)
    {
        _owned.insert(_owned.end(), owned[index] -> begin(), owned[index] -> end());
    }
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: setMinMax                                                                        *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
void DatasetDist::setMinMax(void)
{
    if (_owned.size() == 0) return;

    for(uint index = 0; index < _owned.size(); index += 1)
    {
        updateMinMax(index);
    }
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: setSamples                                                                       *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
void DatasetDist::setSamples(const matrixd &x, const vectord &y)
{
    // Update all data
    _alldata -> setSamples(x,y);

    // Update Owned data
    _owned.clear();

    for(uint index = 0; index < x.size1(); index += 1)
    {
        _owned.push_back(index);
    }

    // Assert verification
    if (_alldata -> getNSamples() != getNSamples()) exit(-1);

    // Update Min and Max indexes
    for(uint index = 0; index < _owned.size(); index += 1)
    {
        updateMinMax(index);
    }
};


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: addSample                                                                        *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
void DatasetDist::addSample(const vectord &x, double y)
{
    // Update all data
    _alldata -> addSample(x,y);

    // Update owned data
    _owned.push_back( (_alldata -> getNSamples() - 1) );

    // Update maximum and minimum values
    updateMinMax(getNSamples() - 1);
}





/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: updateMinMax                                                                     *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
void DatasetDist::updateMinMax(size_t i)
{
    if      ( getSampleY(_min) > getSampleY(i) ) { _min = i; }
    else if ( getSampleY(_max) < getSampleY(i) ) { _max = i; }
};


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: updateMinMax                                                                     *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
void DatasetDist::plotData(TLogLevel level)
{
    // For logging purpose
    FILE_LOG(level) << "Initial points:" ;

    for (uint index = 0; index < getNSamples(); index += 1)
    {
	   FILE_LOG(level) << "X:"  << getSampleX(index)
			           << "|Y:" << getSampleY(index);
    }

    const double  yPoint = getValueAtMinimum();
    const vectord xPoint = getPointAtMinimum();

    FILE_LOG(level) << "Best point so far:" ;
    FILE_LOG(level) << "X:"  << xPoint
		            << "|Y:" << yPoint;
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: eraseData                                                                        *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
void DatasetDist::eraseData(void)
{
    _owned.clear();
}

/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: getSamplesY                                                                      *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
vectord DatasetDist::getSamplesY(void) const
{
    vectord y;

    for (uint index = 0; index < _owned.size(); index += 1)
    {
        utils::append(y, getSampleY(index));
    }

    return y;
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: getSamplesX                                                                      *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
vecOfvec DatasetDist::getSamplesX(void) const
{
    vecOfvec x;

    for (uint index = 0; index < _owned.size(); index += 1)
    {
        x.push_back(getSampleX(index));
    }

    return x;
}


/**************************************************************************************************
 *  Procecure                                                                                     *
 *                                                                                                *
 *  Description: calculateUncertainty                                                              *
 *  Class      : DatasetDist                                                                      *
 **************************************************************************************************/
double DatasetDist::calculateUncertainty(std::vector<uint>& data)
{
    double avg    = 0.0;
    double result = 0.0;

    // Get average
    for (uint index = 0; index < data.size(); index += 1)
    {
        avg += _alldata -> getSampleY(data[index]);
    }

    avg /= data.size();

    // Calculate Uncertainty
    for (uint index = 0; index < data.size(); index += 1)
    {
        result += pow( avg - _alldata -> getSampleY(data[index]), 2);
    }

    return result;
}


} // END of namespace bayesopt
