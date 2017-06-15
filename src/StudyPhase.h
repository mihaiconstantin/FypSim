#ifndef STUDYPHASE_H
#define STUDYPHASE_H


class StudyPhase
{

private:
    Rcpp::NumericMatrix thetaLevels;
    Rcpp::NumericMatrix cellResults;


public:
    StudyPhase(const double &low_level, const double &high_level, const double &increments, const unsigned int &iterations);

    void BuildThetaLevelsMatrix(const Rcpp::DoubleVector &thetaLevelsVector, const unsigned int &replications);

    void ApplyLz(const Rcpp::NumericMatrix &shifted_or_population_parameters, const Rcpp::NumericMatrix &estimated_parameters);

    Rcpp::NumericVector AggregatedIndicators(const Rcpp::NumericMatrix &lz_stats, const Rcpp::NumericVector &thetaLevel, const Rcpp::NumericVector &wlmThetaLevel);

    

    // region getters

    Rcpp::NumericMatrix getThetaLevels() const;

    Rcpp::NumericMatrix getCellResults() const;

    // endregion

};


#endif
