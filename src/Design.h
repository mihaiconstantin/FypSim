#ifndef DESIGN_H
#define DESIGN_H

class Design
{

private:
    Rcpp::DoubleVector shiftProportions;
    Rcpp::DoubleVector shiftMagnitudes;
    Rcpp::NumericVector shiftTypes;
    Rcpp::NumericVector parametersTypes;
    Rcpp::NumericVector testLengths;


public:
    Design(Rcpp::DoubleVector shift_proportions, Rcpp::DoubleVector shift_magnitudes, Rcpp::NumericVector shift_types, Rcpp::NumericVector parameters_types, Rcpp::NumericVector test_lengths);

	Rcpp::NumericMatrix buildDesign();

};

#endif