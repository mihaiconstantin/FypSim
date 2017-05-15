#include <Rcpp.h>
#include "EnvironmentCallsR.h"



// Setting the seed within the R environment. This means
// that the all the seeds the user sets in his R code
// will be the same with this one we use internally
// in our the C++.
void EnvironmentCallsR::SetSeed(const unsigned int &seed)
{
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed = base_env["set.seed"];
    set_seed(seed);
}



// Call R's sequence function and return a double vector we can handle in C++.
// Useful to dynamically determine the levels of theta.
Rcpp::DoubleVector EnvironmentCallsR::Sequence(const double &low_level, const double &high_level, const double &increments)
{
    Rcpp::Environment base_env("package:base");
    Rcpp::Function sequence = base_env["seq"];
    return sequence(low_level, high_level, increments);
}