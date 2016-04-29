/**************************************************************************************************
 *  File:    datasetdist.h                                                                        *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/

#ifndef __DATASETDIST_H__
#define __DATASETDIST_H__

/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include <vector>            // Global

#include "log.hpp"
#include "specialtypes.hpp"  // Bayesopt
#include "dataset.hpp"


namespace bayesopt
{

/**************************************************************************************************
 *  Class: DatasetDist                                                                            *
 **************************************************************************************************/
class DatasetDist : public Dataset
{
public:
    // Constructors
    DatasetDist(Dataset*     data);
    DatasetDist(DatasetDist* data, std::vector<uint>&              owned);
    DatasetDist(DatasetDist* data, std::vector<std::vector<uint>*> owned);

    // Destructors
    virtual ~DatasetDist(){};

    // Methods (redifined)
    void               setSamples       (const matrixd& x    , const vectord& y);
    void               addSample        (const vectord& x    , double y        );
    double             getSampleY       (size_t         index                  ) const;
    vectord            getSampleX       (size_t         index                  ) const;
    double             getLastSampleY   (void                                  ) const;
    vectord            getLastSampleX   (void                                  ) const;
    vectord            getPointAtMinimum(void                                  ) const;
    double             getValueAtMinimum(void                                  ) const;
    size_t             getNSamples      (void                                  ) const;
    void               updateMinMax     (size_t         i                      )      ;
    void               plotData         (TLogLevel      level                  )      ;

    // Methods (new)
    void               eraseData           (void                    )      ;
    std::vector<uint>  getOwnedData        (void                    )      ;
    std::vector<uint>* getOwnedDataPtr     (void                    )      ;
    uint               getOwnedData        (uint               index)      ;
    vectord            getSamplesY         (void                    ) const;
    vecOfvec           getSamplesX         (void                    ) const;
    double             calculateUncertainty(void                    )      ;
    double             calculateUncertainty(std::vector<uint>& data )      ;
    void               setMinMax           (void                    )      ;

protected:
    // Atributes
    Dataset*           _alldata;
    std::vector<uint>  _owned;

    uint               _min;
    uint               _max;
};

/**************************************************************************************************
 *  Inline methods : DatasetDist                                                                  *
 **************************************************************************************************/
inline double             DatasetDist::getSampleY          (size_t index) const { return   _alldata -> mY(_owned[index])  ; }
inline vectord            DatasetDist::getSampleX          (size_t index) const { return   _alldata -> mX[_owned[index]]  ; }
inline double             DatasetDist::getLastSampleY      (void        ) const { return   _alldata -> mY(_owned.back())  ; }
inline vectord            DatasetDist::getLastSampleX      (void        ) const { return   _alldata -> mX[_owned.back()]  ; }
inline vectord            DatasetDist::getPointAtMinimum   (void        ) const { return   _alldata -> getPointAtMinimum(); }
inline double             DatasetDist::getValueAtMinimum   (void        ) const { return   _alldata -> getValueAtMinimum(); }
inline size_t             DatasetDist::getNSamples         (void        ) const { return   _owned.size()                  ; }
inline std::vector<uint>  DatasetDist::getOwnedData        (void        )       { return   _owned                         ; }
inline std::vector<uint>* DatasetDist::getOwnedDataPtr     (void        )       { return (&_owned)                        ; }
inline uint               DatasetDist::getOwnedData        (uint   index)       { return   _owned[index]                  ; }
inline double             DatasetDist::calculateUncertainty(void        )       { return calculateUncertainty(_owned)     ; }

} //namespace bayesopt

#endif // _DATASETDIST_H_
