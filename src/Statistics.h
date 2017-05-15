#ifndef STATISTICS_H
#define STATISTICS_H

class Statistics
{

public:
    static Rcpp::NumericMatrix SimulateGrm(const Rcpp::DoubleVector &theta, const Rcpp::NumericMatrix &item_parameters);

    static Rcpp::NumericMatrix EstimateParametersGrm(const Rcpp::NumericMatrix &data);

    static Rcpp::NumericVector EstimateThetaWml(const Rcpp::NumericVector &data, const Rcpp::NumericMatrix &ipar, double thStart = 0, double tolerance = 1e-9);

    static Rcpp::NumericVector EstimatePolyLzStar(const Rcpp::NumericVector &data, const Rcpp::NumericMatrix &ipar, const double &theta);

    static Rcpp::NumericVector EstimateThetaBisection(Rcpp::NumericVector data, Rcpp::NumericMatrix ipar, double estmeth = 1, double thStart = 0, double tolerance = 1e-5, double bislowin = -4, double bisupin = 4, int nbisin = 100);


    // Wrappers to apply the response-vector-specific functions over an entire dataset.

    static Rcpp::NumericMatrix FullThetaWml(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &item_parameters);

    static Rcpp::NumericMatrix FullPolyLzStar(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &item_parameters, const Rcpp::DoubleVector &theta);

};


#endif
