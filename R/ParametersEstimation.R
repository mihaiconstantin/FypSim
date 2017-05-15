
#' Estimate the item paramters
#' 
#' Fit the GRM using \code{mirt} behind the scenes and 
#' return the simplified parameters. The output will 
#' hold the slope as the first column and from the
#' second the item steps.
#' 
#' @param data (numeric matrix) The data to fit the GRM on.
#' 
#' @export
#' 
estimateParametersGrm <- function(data)
{
	mirt_object = mirt(as.data.frame(data), 1, verbose = FALSE)
	return(coef(mirt_object, IRTpars = TRUE, simplify = TRUE)$items)
}