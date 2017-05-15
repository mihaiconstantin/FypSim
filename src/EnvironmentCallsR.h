#ifndef ENVIRONMENTCALLSR_H
#define ENVIRONMENTCALLSR_H


class EnvironmentCallsR
{

public:
    static void SetSeed(const unsigned int &seed);

    static Rcpp::DoubleVector Sequence(const double &low_level, const double &high_level, const double &increments);
};


#endif