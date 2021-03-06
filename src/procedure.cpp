#include <Rcpp.h>
#include "Design.h"
#include "DesignProcedure.h"
#include "Statistics.h"



//' Builds the factorial design
//' 
//' Builds the factorial design based on the vectors of independent variables.
//' 
//' @param shift_proportions (numeric vector) Proportions of response shift.
//' @param shift_magnitudes (numeric vector) Magnitudes of response shift.
//' @param shift_types (numeric vector) Types of response shift (i.e., towards slope, threshold, or both).
//' @param parameters_types (numeric vector) Types of item parameters (i.e., homogeneous or heterogeneous).
//' @param test_lengths (numeric vector) Lengths of the test (i.e., number of items).
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix buildDesign(Rcpp::DoubleVector shift_proportions,
                                Rcpp::DoubleVector shift_magnitudes,
                                Rcpp::NumericVector shift_types,
                                Rcpp::NumericVector parameters_types,
                                Rcpp::NumericVector test_lengths)
{
    Design Design(shift_proportions, shift_magnitudes, shift_types, parameters_types, test_lengths);
    return Design.buildDesign();
}



//' Applies procedure for a single cell
//' 
//' Performs the calibration and estimation for a single cell. Subsequently, we will use this
//' when we loop through the entire factorial design (i.e., see the runDesign() function).
//' We have to also keep an eye on the fact that the entire design in itself will also
//' need to be replicated \code{n} times.
//' 
//' @param shift_proportion (double) Proportion of response shift.
//' @param shift_magnitude (double) Magnitude of response shift.
//' @param shift_type (int) Type of response shift (i.e., towards slope, threshold, or both).
//' @param parameters_type (int) Type of item parameters (i.e., homogeneous or heterogeneous).
//' @param test_length (int) Lengths of the test (i.e., number of items).
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix runCell(double shift_proportion,
                            double shift_magnitude,
                            unsigned int shift_type,
                            unsigned int parameters_type,
                            unsigned int test_length)
{
    DesignProcedure runProcedure(-2, 2, .5, 9, 500, 500);
    return runProcedure.RunCell(shift_proportion, shift_magnitude, shift_type, parameters_type, test_length);
}



//' Applies procedure for a selected number of cells
//'
//' Performs the calibration and estimation for a selected number of cells.
//' The cells are passed as a numeric matrix where the columns are in the
//' following order: shift proportion, shift magnitude, shift type, the
//' type of the item parameters and the test length. Note, the cells
//' are ran only once. For a design replication we need loop over
//' this function. Thus, we'll have a separate wrapper for that.
//' As a return type, we have the following:
//' Rcpp list containing three numeric matrices: detection rate,
//' mean, and standard deviation. Each of these matrices are
//' structure the following way: the columns indicate the
//' theta level and row numbers points to the same row
//' from withing the selected cells matrix. Thus,
//' values in row \code{n} in any of these
//' matrices represents the value for
//' configuration at row \code{n}
//' in the selected cells.
//'
//' @param selected_cells (numeric matrix) The selected cell configurations to be ran, where
//' the columns respect the order indicated above.
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List runSelectedCells(Rcpp::NumericMatrix selected_cells)
{
    DesignProcedure runProcedure(-2, 2, .5, 9, 500, 500);

    return runProcedure.RunSelectedCells(selected_cells);
}



//' Applies and replicates the procedure for a selected number of cells
//'
//' Very similar with the runSelectedCells() function, the only difference
//' is that it also replicate the selected cells \code{n} times.
//'
//' @param selected_cells (numeric matrix) The selected cell configurations to be ran, where
//' the columns respect the order indicated above.
//' @param design_replications (int) The number of times the selected cells will be replicated.
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List runSelectedCellsWithReplication(Rcpp::NumericMatrix selected_cells, unsigned int design_replications)
{
    DesignProcedure runProcedure(-2, 2, .5, 9, 500, 500);
    return runProcedure.RunSelectedCellsWithReplication(selected_cells, design_replications);
}



//' Compute Cohen's effect size
//'
//' Compute the \code{d} effect size measure for two numeric vectors.
//'
//' @param firstVector (numeric vector) The first vector used in the effect size computations.
//' @param secondVector (numeric vector) The second vector used in the effect size computations.
//'
//' @export
//'
// [[Rcpp::export]]
double effectSizeCohen(Rcpp::NumericVector &firstVector, Rcpp::NumericVector &secondVector)
{
    return Statistics::CohenEffectSize(firstVector, secondVector);
}
