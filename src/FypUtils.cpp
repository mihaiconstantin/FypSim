#include <Rcpp.h>
#include "FypUtils.h"


// Check if a single value passed is either NA or NaN.
bool FypUtils::IsMissing(const double &value)
{
    bool isMissing = false;

    if (R_IsNaN(value) || R_IsNA(value) || ISNA(value) || ISNAN(value) || std::isnan(value))
    {
        isMissing = true;
    }

    return isMissing;
}


// Check if a vector passed contains either NA or NaN.
bool FypUtils::IsMissing(const Rcpp::NumericVector &vectorValues)
{
    for (int i = 0; i < vectorValues.length(); ++i)
    {
        if (Utils::IsMissing(vectorValues(i)))
        {
            return true;
        }
    }

    return false;
}


// Check if a vector passed contains either NA or NaN.
bool FypUtils::IsMissing(const Rcpp::NumericMatrix &matrixValues)
{
    for (int column = 0; column < matrixValues.ncol(); ++column)
    {
        for (int row = 0; row < matrixValues.nrow(); ++row)
        {
            if (Utils::IsMissing(matrixValues(row, column)))
            {
                return true;
            }
        }
    }

    return false;
}
