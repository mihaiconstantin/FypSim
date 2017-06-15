#include <Rcpp.h>
#include "DesignProcedure.h"
#include "Cell.h"
#include "Calibration.h"
#include "StudyPhase.h"
#include <iostream>
#include <time.h>



// Class constructor performing action that are generally applicable to all the
// runner methods that will be called on the instantiated object. For example,
// we can use it to specify the parameters under the cells will be ran.
DesignProcedure::DesignProcedure(const double &start_level,
                                 const double &end_level,
                                 const double &level_increments,
                                 const unsigned int &total_levels,
                                 const unsigned int &level_replications,
                                 const unsigned int &population_sample_size)
{
    // Set the configuration for the cells.
    startLevel = start_level;
    endLevel = end_level;
    levelIncrements = level_increments;
    totalLevels = total_levels;
    levelReplications = level_replications;
    populationSampleSize = population_sample_size;


    // Other generally applicable logic here.


    // region feedback
    // General setup feedback.
    std::cout << "\nFor "
              << totalLevels
              << " levels of theta ("
              << start_level <<  " to "
              << end_level
              << " by "
              << level_increments
              << ") with "
              << levelReplications
              << " replications at each level."
              << std::endl;

    std::cout << "Population sample size: "
              << populationSampleSize
              << "."
              << std::endl;
    // endregion

}



// Runs the study procedure for a single cell and returns a Rcpp
// numeric matrix containing on the rows the detection rate,
// mean, and standard deviation. The columns indicate the
// theta level for which the procedure was applied.
Rcpp::NumericMatrix DesignProcedure::RunCell(const double &shift_proportion,
                                       const double &shift_magnitude,
                                       const unsigned int &shift_type,
                                       const unsigned int &parameters_type,
                                       const unsigned int &test_length)
{
    // Specify the configuration of the cell that is about to be ran.
    Cell Configuration;
    Configuration.shiftProportion = shift_proportion;
    Configuration.shiftMagnitude = shift_magnitude;
    Configuration.shiftType = shift_type;
    Configuration.parametersType = parameters_type;
    Configuration.testLength = test_length;


    // Perform the calibration using the population theta and the current configuration.
    // But first draw a new sample of population theta (for each cell ran).
    const Rcpp::DoubleVector population_theta = Rcpp::rnorm(populationSampleSize, 0, 1);
    Calibration calibratedCell(population_theta, &Configuration);


    // Perform the core of the study (i.e., apply LZ, but make sure you choose the parameters right).
    StudyPhase studiedCell(startLevel, endLevel, levelIncrements, levelReplications);


    // Determine under which model the cell falls (i.e., NULL or ALTERNATIVE) and pick the right item parameters.
    if ((int) shift_proportion > 0)
    {
        studiedCell.ApplyLz(calibratedCell.getShiftedParameters(), calibratedCell.getEstimatedParameters());
    }
    else
    {
        studiedCell.ApplyLz(calibratedCell.getPopulationParameters(), calibratedCell.getEstimatedParameters());
    }


    // Return the results for the current cell.
    return studiedCell.getCellResults();
}



// Runs the study procedure for a selected number cell and returns a
// Rcpp list containing three numeric matrices: detection rate,
// mean, and standard deviation. Each of these matrices are
// structure the following way: the columns indicate the
// theta level and row numbers points to the same row
// from withing the selected cells matrix. Thus,
// values in row \code{n} in any of these
// matrices represents the value for
// configuration at row \code{n}
// in the selected cells.
Rcpp::List DesignProcedure::RunSelectedCells(const Rcpp::NumericMatrix &selected_cells)
{
    // Prepare the "aggregated" matrices.
    Rcpp::NumericMatrix detections(selected_cells.nrow(), totalLevels);
    Rcpp::NumericMatrix means(selected_cells.nrow(), totalLevels);
    Rcpp::NumericMatrix sds(selected_cells.nrow(), totalLevels);
    Rcpp::NumericMatrix effect(selected_cells.nrow(), totalLevels);


    // Get the results for each cell in the selected cells matrix
    // and store the results in the "aggregated" matrices.
    for (int cell = 0; cell < selected_cells.nrow(); ++cell)
    {
        // region feedback
        // Print the current cell number and configuration.
        std::cout << "\tCell "
                  << cell + 1
                  << ": ("
                  << "shift proportion: "   << selected_cells(cell, 0)
                  << " | shift magnitude: " << selected_cells(cell, 1)
                  << " | shift type: "      << selected_cells(cell, 2)
                  << " | parameters type: " << selected_cells(cell, 3)
                  << " | test length: "     << selected_cells(cell, 4)
                  << ") ";
        // endregion

        Rcpp::NumericMatrix cell_data = RunCell(selected_cells(cell, 0), selected_cells(cell, 1), selected_cells(cell, 2), selected_cells(cell, 3), selected_cells(cell, 4));

        detections(cell, Rcpp::_) = cell_data(0, Rcpp::_);
        means(cell, Rcpp::_) = cell_data(1, Rcpp::_);
        sds(cell, Rcpp::_) = cell_data(2, Rcpp::_);
        effect(cell, Rcpp::_) = cell_data(3, Rcpp::_);


        // region feedback
        // Marking the cell competition.
        std::cout << " |done|" << std::endl;
        // endregion
    }

    return Rcpp::List::create(detections, means, sds, effect);
}



// Very similar with the RunSelectedCells() function, the only noticeable difference being that
// this current function does also replicate the selected cells N number of times.
Rcpp::List DesignProcedure::RunSelectedCellsWithReplication(const Rcpp::NumericMatrix &selected_cells, const unsigned int &design_replications)
{
    // region feedback
    // Recording the starting time of the design replications.
    time_t starting_time;
    time(&starting_time);

    std::cout << "Design replications requested: "
              << design_replications
              << " (for "
              << selected_cells.nrow()
              << " cells)."
              << std::endl;

    std::cout << "\n-----------------------------\n"
              << "- Starting the replications."
              << "\n-----------------------------\n"
              << std::endl;
    // endregion


    Rcpp::List replications(design_replications);

    for (unsigned int replication = 0; replication < design_replications; ++replication)
    {
        // region feedback
        // Prints the current replication.
        std::cout << "\nReplication: "
                  << replication + 1
                  << std::endl;
        // endregion

        replications(replication) = RunSelectedCells(selected_cells);
    }


    // region feedback
    // Recording the ending time of the design replications.
    time_t ending_time;

    std::cout << "\n\n-----------------------------------------------\n"
              << "- Time summary: from "
              << starting_time
              << " to "
              <<  time(&ending_time)
              << "."
              << "\n-----------------------------------------------\n"
              << std::endl;
    // endregion

    return replications;
}

