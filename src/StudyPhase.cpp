#include <Rcpp.h>
#include "StudyPhase.h"
#include "EnvironmentCallsR.h"
#include "Statistics.h"
#include "FypUtils.h"
#include <iostream>


// Each cell will perform for a series of theta levels.
// Also, for each level there will be n number of
// performances. We set this configuration now.
StudyPhase::StudyPhase(const double &low_level, const double &high_level, const double &increments, const unsigned int &iterations)
{
    // Prepare the matrix containing the theta levels.
    BuildThetaLevelsMatrix(EnvironmentCallsR::Sequence(low_level, high_level, increments), iterations);
}



// Prepare the matrix where we will store the cell results.
// It goes by the following format: first row holds the
// detection rate, second holds the mean, third holds
// the standard deviation. The columns represent
// the levels of theta selected. Moving next.
//
// Apply the LZ over the matrix containing the levels of theta and
// store the result in the currently instantiated object. Note
// that the parameters involved in this procedure differ on
// the model we run the procedure for (i.e., NULL model
// vs. ALTERNATIVE model). In order to figure out
// what model the cell runs for we can use the
// proportion of shift as an indicator. If
// the proportion of items that will be
// injected is 0, then we know for
// sure we are in the NULL model.
// Otherwise, we are in the
// ALTERNATIVE model.
//
// Now that we established how to differentiate between the NULL and the
// ALTERNATIVE model, it is the time to discuss how core of the paper
// (i.e., the application of the Lz polychotomous star) is used for
// each of the above mentioned scenarios we are interested in.
//
// If we are in the NULL model, FOR EACH REPLICATED LEVEL, we simulate the data using
// the population parameters, estimate the theta via WML method using the estimated
// parameters, and we compute the Lz polychotomous star using the estimated set
// of parameters and the theta WML.
//
// If we are in the ALTERNATIVE model, FOR EACH REPLICATED LEVEL, we simulate the data
// using the shifted parameters, estimate the theta WML using the estimated set of
// parameters, and we compute the Lz polychotomous start using the estimated set
// of parameters and the theta WML.
//
// It actually makes sense. Take a step back and think a bit about it. What we are doing
// is to inject response shift on the data we simulate based on the replications of
// the levels of the theta (we use this because estimating the LZ for just once
// for, say value 3, isn't enough; we want to see how estimating LZ for value
// 3 for 500 times looks like--by taking aggregated statistics for that
// particular level). As a result of the shift injection, we end up
// with a biased dataset. So, what is left to do is to use the
// estimated parameters from the calibration phase (because
// in a real-life scenario we never know the population
// parameters) to simulate some thetas by WML that we
// later use to compute the LZ. Because we computed
// the LZ using the estimated parameters, but we
// apply it on the biased dataset built with the
// shifted parameters, what the LZ flags as
// misfitting responses is the biased we
// introduced! This is, of course,
// valid for the ALTERNATIVE
// model.
void StudyPhase::ApplyLz(const Rcpp::NumericMatrix &shifted_or_population_parameters, const Rcpp::NumericMatrix &estimated_parameters)
{
    cellResults = Rcpp::NumericMatrix(3, thetaLevels.ncol());

    for (int level = 0; level < thetaLevels.ncol(); ++level)
    {
        // The data is estimated either using the shifted parameters (ALTERNATIVE model)
        // or using the population parameters (NULL model).
        Rcpp::NumericMatrix data = Statistics::SimulateGrm(thetaLevels(Rcpp::_, level), shifted_or_population_parameters);


        Rcpp::NumericMatrix theta = Statistics::FullThetaWml(data, estimated_parameters);


        Rcpp::NumericMatrix lz = Statistics::FullPolyLzStar(data, estimated_parameters, theta(Rcpp::_, 0));


        cellResults(Rcpp::_, level) = AggregatedIndicators(lz, thetaLevels(Rcpp::_, level), theta);


        // For safety reasons, let's check every single vector/ matrix we compute to see if NA/ NaN found their way in.
        if(FypUtils::IsMissing(data))        std::cout << " >>> Error: NA | NaN (file: StudyPhase.cpp). Missing in 'data' under 'ApplyLz method'. <<< ";
        if(FypUtils::IsMissing(theta))       std::cout << " >>> Error: NA | NaN (file: StudyPhase.cpp). Missing in 'theta' under 'ApplyLz method'. <<< ";
        if(FypUtils::IsMissing(lz))          std::cout << " >>> Error: NA | NaN (file: StudyPhase.cpp). Missing in 'lz' under 'ApplyLz method'. <<< ";
        if(FypUtils::IsMissing(cellResults)) std::cout << " >>> Error: NA | NaN (file: StudyPhase.cpp). Missing in 'cellResults' under 'ApplyLz method'. <<< ";

    }
}



// Build the matrix containing the theta levels and their appropriate
// number of replications. The rows reflect the number of times the
// theta level is replicated, and the columns represent the
// actual values of the selected theta levels.
void StudyPhase::BuildThetaLevelsMatrix(const Rcpp::DoubleVector &thetaLevelsVector, const unsigned int &replications)
{
    thetaLevels = Rcpp::NumericMatrix(replications, (int) thetaLevelsVector.length());

    for (int level = 0; level < thetaLevelsVector.length(); ++level)
    {
        thetaLevels(Rcpp::_, level) = Rcpp::NumericVector(replications, thetaLevelsVector(level));
    }
}



// Extract the relevant stats we are interested in from the numeric matrix
// resulted aster we estimated the lz. Also, compute the effect size.
Rcpp::NumericVector StudyPhase::AggregatedIndicators(const Rcpp::NumericMatrix &lz_stats, const Rcpp::NumericVector &thetaLevel, const Rcpp::NumericVector &wlmThetaLevel)
{
    Rcpp::LogicalVector detections = lz_stats(Rcpp::_, 1) <= -1.645;

    double detection = Rcpp::mean(detections);
    double mean = Rcpp::mean(lz_stats(Rcpp::_, 1));
    double sd = Rcpp::sd(lz_stats(Rcpp::_, 1));
    double effect = Statistics::CohenEffectSize(thetaLevel, wlmThetaLevel);

    return Rcpp::NumericVector::create(detection, mean, sd, effect);
}



// region getters

Rcpp::NumericMatrix StudyPhase::getThetaLevels() const
{
    return thetaLevels;
}


Rcpp::NumericMatrix StudyPhase::getCellResults() const
{
    return cellResults;
}



// endregion