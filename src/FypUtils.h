#ifndef FYPUTILS_H
#define FYPUTILS_H


class FypUtils
{

public:
    static bool IsMissing(const double &value);
    static bool IsMissing(const Rcpp::NumericVector &vectorValues);
    static bool IsMissing(const Rcpp::NumericMatrix &matrixValues);
};


#endif
