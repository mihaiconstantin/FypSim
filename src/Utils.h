#ifndef UTILS_H
#define UTILS_H


class Utils
{

public:
    static bool IsMissing(const double &value);
    static bool IsMissing(const Rcpp::NumericVector &vectorValues);
    static bool IsMissing(const Rcpp::NumericMatrix &matrixValues);
};


#endif
