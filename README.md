# FypSim
Simulation for the First-Year Paper at Tilburg University (2016-2017) - Person fit statistics

***Work in progress...***

## Sample code

```r
# Specifing the factors.
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


View(design)

# Running a single cell.
# single = runCell(shift_proportion = 50, 
# 		  		  shift_magnitude = .75, 
# 		      		   shift_type = 2,
# 		 		  parameters_type = 0,
# 		      		  test_length = 10)


# Running 3 cells.
selected = runSelectedCells(design[1:10, ])


# Printing results
# single

design[1:10, ]
selected

```
