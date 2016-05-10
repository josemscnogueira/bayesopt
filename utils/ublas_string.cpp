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

        std::string ublas_toString(const matrixd& mat)
        {
            std::stringstream output;
                              output << mat;

            return output.str();
        }

        vectord string_toUblasVectord(const std::string& s)
        {
            std::stringstream input;
                              input << s;

            vectord result;

            input >> result;

            return result;
        }

        matrixd string_toUblasMatrixd(const std::string& s)
        {
            std::stringstream input;
                              input << s;

            matrixd result;

            input >> result;

            return result;
        }
    } // End of namespace utils
} // End of namespace bayesopt
