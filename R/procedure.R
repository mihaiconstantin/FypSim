# 
# # # # # # # # # # # # # # # # # # # # #
# # Script for the FYP data (simulations)
# # # # # # # # # # # # # # # # # # # # #
# 
# 
# # We start by loading the necessary packages.
# 
# 
# library(Rcpp)
# library(mirt)
# library(changeIRT)
# source("D:/01_School/01_Faculty/01_My Psy/Master/Tilburg University/00_Princeps/_archive/1st Year/Internal Traineeship I/_assignments/_R code/package/ChangeIRT/R/temp/PolyLzStar.R")
# source("D:/01_School/01_Faculty/01_My Psy/Master/Tilburg University/00_Princeps/_archive/1st Year/Internal Traineeship I/_assignments/_R code/package/ChangeIRT/R/temp/FsimGRM.R")
# 
# 
# # We have a factorial design, so it might be handy to
# # define all our conditions beforehand, then code
# # them accordingly. A solution is to use store
# # all the conditions in vectors and keep in
# # mind the what represents what. This way
# # we can build a matrix and handle the
# # simulations for the paper easier.
# 
# 
# # Let's build our factorial design and describe
# # the conditions we have, so they make sense
# # for everybody who wants to read the code.
# 
# # Test lenght (3 levels):
# # 	- 5 items
# #	- 10 items
# #	- 30 items
# 
# # Item parameters (2 levels):
# # 	- 0 homogenous items
# # 	- 1 heterogenous items
# 
# # Response shift proportion of the test lenght (4 levels):
# # 	- 0%
# # 	- 10%
# #	- 20%
# #	- 30%
# #
# 
# # Response shift magnitude towards the item parameters (4 levels)
# #	- -.75
# #	- -.25
# #	- .25
# #	- .75
# 
# # Response shift direction (3 levels):
# # 	- 0 (response shift bias towards slope)
# # 	- 1 (response shift bias towards threshold)
# #	- 2 (response shift bias towards both the slope and threshold)
# 
# 
# # Bases on the following combiation of factors and the research
# # question, we can derive two conditions of interest for the
# # study. First, when the response shift propotion is 0 for
# # all the other combinations are are under the control
# # condition (NULL model). Second, when the response
# # shift is different than 0, then we are under
# # the experimental condition, which we call
# # the ALTERNATIVE model.
# 
# 
# # Let's discuss few words about the procedure and the indicators we
# # want to analize. Under the NULL model we want to look at the
# # Type I error rates. Under the ALTERNATIVE model we want to
# # look at the detection rates (power). For each of these,
# # we will look at the mean, standard deviation, and the
# # skeewness and kurtosis indicators. This is all done
# # at an agragated level, meaning that the entire
# # design is replicated 100 times, and we do
# # not take the mean of just one single
# # run (aka, simulation of data).
# 
# 
# # We want to be able to use this in a really easy way, that
# # is, we want to have a vector for each of the factors we
# # enumarated. Then, we want to be able to throw them in
# # a faction and get the indicators for both models.
# # Then we will run that function 100 times and
# # get our agregated indicators and build
# # some graphs and presentation tables.
# 
# 
# #
# 
# 
# 
# 
# 
# 
# # Define the factor vectors.
# 
# shift.prop = c(0, 10, 20, 30)
# shift.magn = c(-.75, -.25, .25, .75)
# shift.type = c(0, 1, 2)
# test.lenght = c(5, 10, 30)
# item.para = c(0, 1)
# 
# 
# # Define the values for the item parameters.
# 
# # Homogenous items.
# 
# sample.ho <- function(j)
# {
# 	a.ho = runif(j, 1.5, 2.5)
# 	b.ho = outer(runif(j, 0, 1.25), c(-1, -.5, .5, 1), "+")
# 
# 	return(list(a.ho, b.ho))
# }
# 
# 
# # Heterogeneous items.
# 
# sample.he <- function(j)
# {
# 	a.he = runif(j, 1, 2.5)
# 	b.he = outer(runif(j, -1, 2.5), c(-.5, -.2, .2, .5), "+")
# 
# 	return(list(a.he, b.he))
# }
# 
# 
# a.he = c()
# a.ho = c()
# 
# steps.he = c()
# steps.ho = c()
# 
# 
# 
# 
# # It is also the time to define the levels of theta we want to
# # test the precedre on and determine the amount of decided
# # replications we want to use. In our case we will draw
# # 5 levels of theta from the U(0, 1) and replicate
# # them for 500 times. We store everything in a
# # matrix and use the columns as a vector of
# # replications.
# 
# population_theta = rnorm(500)
# levels = sapply(c(-2, -1, 0, 1, 2), rep, 500)
# 
# 
# # Perform the procedure once with the given factors and
# # return a matrix with two columns and four rows,
# # where representing the indicators for both
# # of models we are interested in (i.e.,
# # NULL vs ALTERNATIVE).
# 
# #run.1 = perform()
# 
# #r
# 
# 
# # done
# 
# 
# 
# 
# 
# 
# # Grup the matrices by null vs alternative
# 
# null = conditions[conditions[, 1] == 0, ]
# alternative = conditions[conditions[, 1] != 0, ]
# 
# null = null[order(null[, 5]), ]
# alternative = alternative[order(alternative[, 5]), ]
# 
# 
# # Testing if sorting is okay.
# 
# null
# 
# alternative
# 
# 
# # Here we start by manipulating the response shift and then we run
# # the procedure proposed. As a result of this, for each i in the
# # loop, which represnts a "cell" in our design, we want to get
# # the indicators for how well lz* perefromed. Later we will
# # think of a nicer way to store these results and compute
# # them in a cleaner way. This is just a first try.
# 
# 
# res.alternative = matrix(NA, 500, dim(alternative)[1])
# 
# 
# # This is the main loop for the theta levels
# 
# for(level in dim(levels)[2]) {
# 
# 	# This is the secondary loop for the item parameters and response shift.
# 
# 	for(i in 1:dim(alternative)[1]) {
# 
# 		# Picking the proper parameters.
# 
# 		if(alternative[i, 4] == 0)
# 		{
# 			set.seed(1993)
# 			param = sample.ho(alternative[i, 5])
# 		} else
# 		{
# 			set.seed(1993)
# 			param = sample.he(alternative[i, 5])
# 		}
# 
# 
# 		# Here we perform the CALIBRATION.
# 
# 		# sim = FsimGRM(population_theta, cbind(param[[1]], param[[2]]))
# 
# 		# sim = simulateirt(population_theta, param[[2]], param[[1]])
# 		sim = FsimGRM(population_theta, cbind(param[[1]], param[[2]]))
# 		est = estimateirt(sim[[2]])
# 
# 
# 		# Now we store the estimated parameters as a result of the CALIBRATION.
# 		param.e = list(est$get.parameters[, 1], est$get.parameters[, 2:5])
# 
# 
# 		# Deciding where and how much shift is injected.
# 		# Short notice, shift will be injected on the
# 		# estimated parameters.
# 
# 		param.s = param
# 		proportion = ceiling(alternative[i, 1] / 100 * alternative[i, 5])
# 
# 		head(alternative, 2)
# 
# 		if(alternative[i, 3] == 0)
# 		{
# 			param.s[[1]][1:proportion] = param[[1]][1:proportion] + alternative[i, 2]
# 
# 		} else if(alternative[i, 3] == 1)
# 		{
# 			param.s[[2]][1:proportion, ] = param[[2]][1:proportion, ] + alternative[i, 2]
# 
# 		} else {
# 			param.s[[1]][1:proportion] = param[[1]][1:proportion] + alternative[i, 2]
# 			param.s[[2]][1:proportion, ] = param[[2]][1:proportion, ] + alternative[i, 2]
# 		}
# 
# 
# 		# Here we perform the STUDY PHASE.
# 
# 		# sim2 = simulateirt(levels[, level], param.s[[2]], param.s[[1]])
# 		sim2 = FsimGRM(levels[, level], cbind(param.s[[1]], param.s[[2]]))
# 
# 
# 		# Let it estimate.
# 		lz = PolyLZstar(sim2[[2]], cbind(param.e[[1]], param.e[[2]]))
# 
# 
# 		# Try to store the results somehow.
# 		res.alternative[, i] = lz[, 3]
# 
# 		if (i == 1) break
# 
# 	} # end of secondary
# 
# 
# 	# Stopping at next theta level.
# 	if (level == 2) break
# 
# } # end of primary
# 
# 
# # Did it work???????
# 
# 
# View(round(res.alternative, 5))
# 
# 
# 
# 
# 
# 
# 
# mean(res.alternative[, 1] < -1.645)
# 
# 
# 
# #
# 
# 
# 
# #