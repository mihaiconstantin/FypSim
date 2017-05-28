#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "Cell.h"
#include "Calibration.h"
#include "Statistics.h"
#include <iostream>



// Initialize the Calibration object and pass the relevant arguments.
// We will also pass the Cell struct by reference.
Calibration::Calibration(const Rcpp::DoubleVector &population_theta, Cell *cell_configuration)
{
	// Store the addresses of the data we'll use for computations.
	cellConfiguration = cell_configuration;


    // Specify the steps intervals. Later we might also want to
    // allow the researcher to pass this as a vector from R.
    steps = Rcpp::DoubleVector::create(-1, -.5, .5, 1);


    // Build the structure for the item parameters.
    populationParameters = Rcpp::NumericMatrix(cellConfiguration->testLength, (unsigned int) steps.length() + 1);
    estimatedParameters = Rcpp::NumericMatrix(cellConfiguration->testLength, (unsigned int) steps.length() + 1);


    // Sample the appropriate population parameters. Then store the estimated parameters
    // based on generated data. Inject the response shift in the right amount, at the
    // right place. Also make sure that the estimated parameters via MIRT do not
    // contain NaNs or NAs. If they do contain, then re-run the calibration
    // phase (i.e., the 3 methods).
    SampleParameters();
    EstimateParameters(population_theta);
    InjectResponseShift();

    bool assumeNaN = true;

    while (assumeNaN)
    {
        for (int column = 0; column < estimatedParameters.ncol(); ++column)
        {
            for (int row = 0; row < estimatedParameters.nrow(); ++row)
            {
                if(std::isnan(estimatedParameters(row, column)))
                {
                    // region feedback
                    std::cout << "\t>>> (!) 'MIRT' estimation produced 'NaNs'. Re-calibrating the cell. <<< ";
                    // endregion

                    SampleParameters();
                    EstimateParameters(population_theta);
                    InjectResponseShift();
                }
                else
                {
                    assumeNaN = false;
                }
            }
        }
    }

}



// Sample either homogeneous or heterogeneous parameters based on
// the value of the parameter_type field for the current cell.
void Calibration::SampleParameters()
{
	// Define the data structures that we will populate later.
	Rcpp::DoubleVector means;


	// In sampling the item parameters, don't forget to seed the seed
	// that the user has passed via the constructor. The reason why
	// we don't need to give SetSeed a value is that we store the
	// the integer passed by the user in a protected field in
	// the class. This way, we are always sure that we use
	// the same seed.
	if (cellConfiguration->parametersType)
	{
        populationParameters(Rcpp::_, 0) = Rcpp::runif(cellConfiguration->testLength, 1.5, 2.5);
		means = Rcpp::runif(cellConfiguration->testLength, 0, 1.25);
	} 
	else 
	{
        populationParameters(Rcpp::_, 0) = Rcpp::runif(cellConfiguration->testLength, 1, 2.5);
		means = Rcpp::runif(cellConfiguration->testLength, -1.5, 2.5);
	}


	// Populate the parameters matrix with the thresholds by
	// iterating over the steps and adding the values as
	// columns to the matrix. Note that the first one
	// has already been set with the a parameters.
	for (int i = 0; i < steps.length(); i++)
	{
        populationParameters(Rcpp::_, i + 1) = means + steps(i);
	}
}



// First simulate data on using the population parameters and the theta.
// Then, account for the sampling error and fit the GRM model on the
// simulated data. Finally, store the estimated parameters and the
// calibration phase is concluded.
void Calibration::EstimateParameters(const Rcpp::DoubleVector &population_theta)
{
    Rcpp::NumericMatrix data = Statistics::SimulateGrm(population_theta, populationParameters);

    estimatedParameters = Statistics::EstimateParametersGrm(data);
}



// Inject the appropriate response shift configuration
// at the item parameters and store a copy of those
// modified values at the item parameters in the
// appropriate class field.
void Calibration::InjectResponseShift()
{
    // Convert the Rcpp numeric matrix to RcppArmadillo for easier cuts.
    arma::mat temp = Rcpp::as<arma::mat>(populationParameters);


    // We must subtract 1 because C++ indices start at 0 (i.e., a regular Rcpp matrix).
    // Therefore, if we get a proportion of 3, meaning 30%, we know that we need to
    // inject shift for the first three items. So, we subtract 1 are we remain
    // with 2. In the context of a matrix of parameters this means rows 0, 1,
    // and 2. Thus, we have successfully inject shift for the first three
    // parameters of interest. See below.
    unsigned int proportion = (unsigned int) ceil(cellConfiguration->shiftProportion / 100 * cellConfiguration->testLength);


    // To avoid scenarios where the user types in wired proportions
    // we safely determine until what item we cut the matrix of
    // parameters to inject sift by checking some conditions.
    if (proportion > 0)
    {
        unsigned int until_item = proportion - 1;
        if (proportion > cellConfiguration->testLength) { until_item = cellConfiguration->testLength - 1; }

        switch (cellConfiguration->shiftType)
        {
            case 0:
                temp.submat(0, 0, until_item, 0) = temp.submat(0, 0, until_item, 0) + cellConfiguration->shiftMagnitude;
                break;

            case 1:
                temp.submat(0, 1, until_item, (unsigned int) steps.length()) = temp.submat(0, 1, until_item, (unsigned int) steps.length()) + cellConfiguration->shiftMagnitude;
                break;

            default:
                temp.submat(0, 0, until_item, 0) = temp.submat(0, 0, until_item , 0) + cellConfiguration->shiftMagnitude;
                temp.submat(0, 1, until_item, (unsigned int) steps.length()) = temp.submat(0, 1, until_item, (unsigned int) steps.length()) + cellConfiguration->shiftMagnitude;
        }
    }

    shiftedParameters = Rcpp::wrap(temp);
}



// region getters


Rcpp::NumericMatrix Calibration::getPopulationParameters() const
{
    return populationParameters;
}


Rcpp::NumericMatrix Calibration::getEstimatedParameters() const
{
    return estimatedParameters;
}


Rcpp::NumericMatrix Calibration::getShiftedParameters() const
{
    return shiftedParameters;
}


Cell *Calibration::getCellConfiguration() const
{
    return cellConfiguration;
}

// endregion