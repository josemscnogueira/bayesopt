#ifndef __UBLAS_STRING_HPP__
#define __UBLAS_STRING_HPP__

#include <string>
#include <sstream>
#include "specialtypes.hpp"

namespace bayesopt
{
    namespace utils
    {
        std::string ublas_toString       (const vectord&     vec);
        std::string ublas_toString       (const matrixd&     mat);
        vectord     string_toUblasVectord(const std::string& s);
        matrixd     string_toUblasMatrixd(const std::string& s);
    } // End of namespace utils
} // End of namespace bayesopt


#endif
