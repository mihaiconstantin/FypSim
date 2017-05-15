#ifndef CALIBRATION_H
#define CALIBRATION_H

class Calibration
{

private:
    struct Cell *cellConfiguration;

    Rcpp::DoubleVector steps;

    Rcpp::NumericMatrix populationParameters;

    Rcpp::NumericMatrix estimatedParameters;

    Rcpp::NumericMatrix shiftedParameters;


public:
	Calibration(Rcpp::DoubleVector &population_theta, Cell *cell_configuration);

	void SampleParameters();

    void EstimateParameters(const Rcpp::DoubleVector &population_theta);

    void InjectResponseShift();

    // region getters

    Rcpp::NumericMatrix getPopulationParameters() const;

    Rcpp::NumericMatrix getEstimatedParameters() const;

    Rcpp::NumericMatrix getShiftedParameters() const;

    Cell *getCellConfiguration() const;

    // endregion
};

#endif