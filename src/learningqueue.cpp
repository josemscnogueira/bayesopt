/**************************************************************************************************
 *  File:    learningqueue.cpp                                                                    *
 *  Author:  Jose Miguel Nogueira, josemscnogueira@gmail.com                                      *
 *                                                                                                *
 *  History:                                                                                      *
 **************************************************************************************************/


/**************************************************************************************************
 *  Include Files                                                                                 *
 **************************************************************************************************/
#include "learningqueue.hpp"
#include <boost/lexical_cast.hpp>

LearningElement::LearningElement(void)
{
    this -> query   = vectord(1,0);
    this -> quality = 0.0;
}


LearningElement::LearningElement(vectord query, double quality)
{
    this -> query         = query;
    this -> quality       = quality;
}


bool LearningElement::operator()(const LearningElement &lhs, const LearningElement &rhs) const
{
    return (lhs.quality < rhs.quality);
}

/**************************************************************************************************
 *                                                                                                *
 *                                                                                                *
 *                                                                                                *
 *                                                                                                *
 *                                                                                                *
 *                                                                                                *
 *                                                                                                *
 *                                                                                                *
 *                                                                                                *
 **************************************************************************************************/

LearningQueueWrapper::LearningQueueWrapper(size_t size)
{
    _size = size;

    clear();
}


void LearningQueueWrapper::clear(void)
{
    while (_queue.size() > 0) _queue.pop();
}


LearningElement LearningQueueWrapper::operator[](uint index)
{
    LearningQueue   aux_queue;
    LearningElement result;

    while (index != 0)
    {
        aux_queue.push(_queue.top());
        _queue   .pop();
        index -= 1;
    }

    result = _queue.top();

    while (aux_queue.size() != 0)
    {
        _queue   .push(aux_queue.top());
        aux_queue.pop();
    }

    return result;
}


void LearningQueueWrapper::resize(uint size)
{
    _size = size;
}


bool LearningQueueWrapper::empty(void)
{
    return _queue.empty();
}


size_t LearningQueueWrapper::size(void)
{
    return _queue.size();
}


LearningElement LearningQueueWrapper::top(void)
{
    return _queue.top();
}


void LearningQueueWrapper::push(LearningElement& elem)
{
    _queue.push(elem);

    while (_queue.size() > _size) _queue.pop();
}


void LearningQueueWrapper::push(vectord query, double quality)
{
    LearningElement new_elem = LearningElement(query, quality);

    push(new_elem);
}


void LearningQueueWrapper::pop(void)
{
    _queue.pop();
}


std::string LearningQueueWrapper::toString(void)
{
    std::string result;
    uint        length = size();

    if (length != 0)
    {
        for (uint index = 0; index < length; index += 1)
        {
            result += boost::lexical_cast<std::string>(-(*this)[length - 1 - index].quality);

            if (index != _size - 1) result += std::string(", ");
        }

        if (length != _size)
        {
            for (uint index = 0; index < (_size - length); index += 1)
            {
                result += boost::lexical_cast<std::string>(0.0);

                if (index != _size - length - 1) result += boost::lexical_cast<std::string>(", ");
            }
        }
    }

    return result;
}


std::vector<vectord> LearningQueueWrapper::getBestQueries(void)
{
    std::vector<vectord> results;

    for (uint index = 0; index < _size; index += 1)
    {
        results.push_back( (*this)[index].query );
    }

    return results;
}