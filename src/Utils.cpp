#include "Utils.h"
#include <R.h>


// Check if the values passed is either NA or NaN.
bool Utils::IsMissing(const double &value)
{
    bool isMissing = false;

    if (R_IsNaN(value) || R_IsNA(value) || ISNA(value) || ISNAN(value))
    {
        isMissing = true;
    }

    return isMissing;
}


// Next method.