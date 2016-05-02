#include "ublas_string.hpp"

namespace bayesopt
{
    namespace utils
    {
        std::string ublas_toString(const vectord& vec)
        {
            std::stringstream output;
                              output << vec;

            return output.str();
        }
    } // End of namespace utils
} // End of namespace bayesopt
