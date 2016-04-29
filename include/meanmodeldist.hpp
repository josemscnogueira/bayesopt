/**************************************************************************************************
 *  File:    meanmodeldist.h                                                                      *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/

#ifndef __MEANMODELDIST_HPP__
#define __MEANMODELDIST_HPP__

/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include <vector>              // Global

#include "specialtypes.hpp"    // Bayesopt
#include "mean_functors.hpp"


namespace bayesopt
{

/**************************************************************************************************
 *  Class: MeanModelDist                                                                          *
 **************************************************************************************************/
class MeanModelDist : public MeanModel
{
public:
    // Constructors
    MeanModelDist(MeanModel*     mean);
    MeanModelDist(MeanModelDist* mean, std::vector<uint>& owned);
    MeanModelDist(MeanModelDist* mean, std::vector<uint>* owned);

    // Destructors
    virtual ~MeanModelDist(){};

    // Methods (redifined)
    ParametricFunction* getMeanFunc  (void                 );
    void                setParameters(const vectord&  theta);
    vectord             getParameters(void                 );
    size_t              nParameters  (void                 );
    vectord             getFeatures  (const vectord&  x    );
    void                getFeatures  (const vectord&  x    , vectord& kx);
    size_t              nFeatures    (void                 );
    void                setPoints    (const vecOfvec& x    );
    void                addNewPoint  (const vectord&  x    );
    vectord             muTimesFeat  (void                 );
    double              muTimesFeat  (const vectord&  x    );
    void                setMean      (const vectord&  muv  , const vectord& smu, std::string m_name, size_t dim);
    void                setMean      (mean_parameters mean , size_t         dim);

    // Methods (new)
    matrixd             getmFeatM    (void);
    uint                getNSamples  (void);



protected:
    // Atributes
    MeanModel*          _alldata;
    std::vector<uint>   _owned;
};

/**************************************************************************************************
 *  Inline methods : MeanModelDist                                                                *
 **************************************************************************************************/
inline ParametricFunction* MeanModelDist::getMeanFunc  (void                              ) { return _alldata -> getMeanFunc    (         ); }
inline void                MeanModelDist::setParameters(const vectord&  theta             ) {        _alldata -> setParameters  (theta    ); }
inline vectord             MeanModelDist::getParameters(void                              ) { return _alldata -> getParameters  (         ); }
inline size_t              MeanModelDist::nParameters  (void                              ) { return _alldata -> nParameters    (         ); }
inline vectord             MeanModelDist::getFeatures  (const vectord&  x                 ) { return _alldata -> getFeatures    (x        ); }
inline void                MeanModelDist::getFeatures  (const vectord&  x    , vectord& kx) {        _alldata -> getFeatures    (x    , kx); }
inline size_t              MeanModelDist::nFeatures    (void                              ) { return _alldata -> nFeatures      (         ); }

inline vectord             MeanModelDist::muTimesFeat  (void                              ) { return boost::numeric::ublas::prod      (_alldata -> getmMu(), getmFeatM()               ); }
inline double              MeanModelDist::muTimesFeat  (const vectord& x                  ) { return boost::numeric::ublas::inner_prod(_alldata -> getmMu(), _alldata -> getFeatures(x)); }

inline void                MeanModelDist::setMean      (const vectord&  muv  , const vectord& smu, std::string m_name, size_t dim) { _alldata -> setMean(muv, smu, m_name, dim); }
inline void                MeanModelDist::setMean      (mean_parameters mean , size_t         dim                                ) { _alldata -> setMean(mean, dim            ); }
inline uint                MeanModelDist::getNSamples  (void                              ) { return _owned.size(); };

} //namespace bayesopt

#endif // _MEANMODELDIST_H_
