#ifndef DESIGNPROCEDURE_H
#define DESIGNPROCEDURE_H


class DesignProcedure
{

private:
    double startLevel;
    double endLevel;
    double levelIncrements;
    unsigned int totalLevels;
    unsigned int levelReplications;
    unsigned int populationSampleSize;

public:
    DesignProcedure(
            const double &start_level,
            const double &end_level,
            const double &level_increments,
            const unsigned int &total_levels,
            const unsigned int &level_replications,
            const unsigned int &population_sample_size);


    Rcpp::NumericMatrix RunCell(const double &shift_proportion,
                                const double &shift_magnitude,
                                const unsigned int &shift_type,
                                const unsigned int &parameters_type,
                                const unsigned int &test_length);


    Rcpp::List RunSelectedCells(const Rcpp::NumericMatrix &selected_cells);


};


#endif
