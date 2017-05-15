#include <Rcpp.h>
#include "Design.h"


// Storing the references to the independent variables.
// We'll use these memory addresses to construct the
// all possible cells for our factorial design.
Design::Design(Rcpp::DoubleVector shift_proportions,
               Rcpp::DoubleVector shift_magnitudes,
               Rcpp::NumericVector shift_types,
               Rcpp::NumericVector parameters_types,
               Rcpp::NumericVector test_lengths)
{
    shiftProportions = shift_proportions;
    shiftMagnitudes = shift_magnitudes;
    shiftTypes = shift_types;
    parametersTypes = parameters_types;
    testLengths = test_lengths;
}



// Using the memory addresses stored in the class
// members, we build the cells of our factorial
// design and output everything as a numeric
// matrix data type (i.e., as in Rcpp).
Rcpp::NumericMatrix Design::buildDesign()
{
	int rows =  (int) shiftProportions.length() * (int) shiftMagnitudes.length() * (int) shiftTypes.length() * (int) parametersTypes.length() * (int) testLengths.length();
	Rcpp::NumericMatrix design(rows, 5);
	int row = 0;

	for (int a = 0; a < shiftProportions.length(); a++)
	{
		for (int b = 0; b < shiftMagnitudes.length(); b++)
		{
			for (int c = 0; c < shiftTypes.length(); c++)
			{
				for (int x = 0; x < parametersTypes.length(); x++)
				{
					for (int y = 0; y < testLengths.length(); y++)
					{
						design(row, Rcpp::_) = Rcpp::NumericVector::create(shiftProportions(a), shiftMagnitudes(b), shiftTypes(c), parametersTypes(x), testLengths(y));
						row++;
					}
				}
			}
		}
	}
	
	return design;
}


// 



