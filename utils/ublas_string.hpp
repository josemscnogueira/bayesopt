#ifndef __UBLAS_STRING_HPP__
#define __UBLAS_STRING_HPP__

#include <string>
#include <sstream>
#include "specialtypes.hpp"

namespace bayesopt
{
    namespace utils
    {
        std::string ublas_toString(const vectord& vec);
        std::string ublas_toString(const matrixd& mat);
    } // End of namespace utils
} // End of namespace bayesopt


#endif
