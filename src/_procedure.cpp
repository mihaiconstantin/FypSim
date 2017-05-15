#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "Design.h"
#include "Calibration.h"
#include "StudyPhase.h"
#include "Cell.h"



//' Build the factorial design
//' 
//' Build the factorial design based on the vectors of independent variables.
//' 
//' @param s_prop (numeric vector) Proportions of response shift.
//' @param s_magn (numeric vector) Magnitudes of response shift.
//' @param s_type (numeric vector) Types of response shift (i.e., towards slope, threshold, or both).
//' @param p_type (numeric vector) Types of item parameters (i.e., homogeneous or heterogeneous).
//' @param t_leng (numeric vector) Lengths of the test (i.e., number of items).
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



//' Apply procedure for a single cell
//' 
//' Perform the calibration and estimation for a single cell. Subsequently, we will use this
//' when we loop through the entire factorial design and for each cell we instantiate the
//' SimulateCell class. We have to also keep an eye on the fact that the entire design
//' in itself will also need to be replicated \code{n} times.
//' 
//' @param population_theta (double) Vector of latent abilities for the population.
//' @param shift_proportion (double) Proportion of response shift.
//' @param shift_magnitude (double) Magnitude of response shift.
//' @param shift_type (int) Type of response shift (i.e., towards slope, threshold, or both).
//' @param parameters_type (int) Type of item parameters (i.e., homogeneous or heterogeneous).
//' @param test_length (int) Lengths of the test (i.e., number of items).
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix runCell(Rcpp::DoubleVector population_theta,
                            double shift_proportion,
                            double shift_magnitude,
                            unsigned int shift_type,
                            unsigned int parameters_type,
                            unsigned int test_length)
{
    // Specify the configuration of the cell that is about to be ran.
    Cell Configuration;
    Configuration.shiftProportion = shift_proportion;
    Configuration.shiftMagnitude = shift_magnitude;
    Configuration.shiftType = shift_type;
    Configuration.parametersType = parameters_type;
    Configuration.testLength = test_length;


    // Perform the calibration using the population theta and the current configuration.
    // TODO: Since we draw new population item parameters each time we run a cell, shouldn't we also draw a new vector of population thetas?
    Calibration calibratedCell(population_theta, &Configuration);


    // Perform the core of the study (i.e., apply LZ, but make sure you choose the parameters right).
    //StudyPhase studiedCell(-2, 2, .5, 50);
    //
    // Determine under which model the cell falls (i.e., NULL or ALTERNATIVE) and pick the right item parameters.
    //if ((int) shift_proportion > 0)
    //{
    //    studiedCell.ApplyLz(calibratedCell.getShiftedParameters(), calibratedCell.getPopulationParameters());
    //}
    //else
    //{
    //    studiedCell.ApplyLz(calibratedCell.getPopulationParameters(), calibratedCell.getPopulationParameters());
    //}
    //
    //
    // Return the results for the current cell.
    //return studiedCell.getCellResults();

    return calibratedCell.getShiftedParameters();

}