/**************************************************************************************************
 *  File:    learningqueue.hpp                                                                    *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/

#ifndef __LEARNINGQUEUE_HPP__
#define __LEARNINGQUEUE_HPP__

#include <vector>            // Global
#include <queue>
#include <string>

#include "specialtypes.hpp"  // Bayesopt


/**************************************************************************************************
 *  Class: LearningElement                                                                        *
 **************************************************************************************************/
class LearningElement
{
public:
    // Attributes
    vectord      query;
    double       quality;

    // Constructor
    LearningElement(void);
    LearningElement(vectord query, double quality);

    // Methods
    bool operator()(const LearningElement &lhs, const LearningElement &rhs) const;
};
typedef std::priority_queue<LearningElement, std::vector<LearningElement>, LearningElement> LearningQueue;


/**************************************************************************************************
 *  Class: LearningQueueWrapper                                                                   *
 **************************************************************************************************/
class LearningQueueWrapper
{
protected:
    // Attributes
    size_t         _size;
    LearningQueue  _queue;

public:
    // Constructor
    LearningQueueWrapper(size_t size = 5);

    // Operator
    LearningElement      operator[](uint index);

    // Methods
    void                 clear         (void);
    void                 resize        (uint             size);
    bool                 empty         (void);
    size_t               size          (void);
    LearningElement      top           (void);
    void                 push          (LearningElement& elem);
    void                 push          (vectord          query, double quality);
    void                 pop           (void);
    std::vector<vectord> getBestQueries(void);
    std::string          toString      (void);
};

#endif