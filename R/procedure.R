# A bunch of functions that allows us to use other people's PCs and simulate
# parts of our factorial design. They call the four main C++ functions
# behind the scenses, but it makes it convenient to save the data
# after each 50 replicated cells.



#' Internal function that initalizes the factorial design from the default values.
#' 
#' @export
#' 
initializeDesign <- function()
{
	shift.prop   = c(0, 10, 20, 50)
	shift.magn   = c(-.75, -.25, .25, .75)
	shift.type   = c(0, 1, 2)
	item.para    = c(0, 1)
	test.length  = c(5, 10, 30)


	# Building the factorial design.
	design = buildDesign(shift.prop,
						 shift.magn,
						 shift.type,
						 item.para,
						 test.length)

	return(design)
}



#' Function that allows to devide the simulations across multiple PCs.
#' It stores both the data for the replications, and the designs.
#' Storing the designs might help us spot possible mistakes.
#' 
#' @export
#' 
runCellRange <- function(start, end, directory = ".", replications = 100)
{
	setwd(directory)
	
	
	designs = initializeDesign()
	results = runSelectedCellsWithReplication(designs[start:end, ], replications)
	
	
	saveRDS(object = results, file = paste(start, "_to_", end, "_results.RData", sep = ""))
	saveRDS(object = designs, file = paste(start, "_to_", end, "_designs.RData", sep = ""))
	
	
	cat(paste(start, "_to_", end, "_results.RData", sep = ""), " was written sucessfuly at ", directory, "\n\n")
}



