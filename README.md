# FypSim

**Convenient data simulations** for a `Person-Fit Statistic` (l<sub>z</sub><sup>*</sup> for polychotomous items) using `R` and `C++`. Part of the First-Year Paper project at Tilburg University (2016-2017) &ndash; Individual Differences and Assessment Research Master's programme.


## Description
The code block below illustrates the main four functions in the package. Using them, one is able to replicate the simulations for the First-Year Paper or test how the procedure performs against a specific combination of factors (i.e., see the package documentation). Also, make sure to check the documentation for each function listed below.

Alongside with the source files (`.cpp` and `.h`), this repository also contains the compiled files (`.dll` and `.o`) for a `64bit` platform. You can also compile for your specific platform. Note that this package relies heavily on the [Rcpp](https://github.com/RcppCore/Rcpp) and [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo) packages. Make sure you have the appropriate tools ready (e.g., see [section 1.3](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf) from the `Rcpp` vignettes). We recommend using the [devtools](https://github.com/hadley/devtools) package for a straightforward installation and compilation, once all the necessary tools are in place.

### Exported functions

#### buildDesign()
- building a matrix based on the combinations of all the factors in the study
```r
# Specifing the factors.
shift.proportion = c(0, 10, 25, 50, 70)
shift.magnitude  = c(-1, -.75, -.25, .25, .75, 1)
shift.type       = c(0, 1, 2)
parameters.type  = c(0, 1)
test.length      = c(5, 10, 30)
 
 
# Building the factorial design.
design = buildDesign(shift.proportion,
		     shift.magnitude,
		     shift.type,
		     parameters.type,
		     test.length)
		     
> head(design)
     [,1]  [,2] [,3] [,4] [,5]
[1,]    0 -0.75    0    0    5
[2,]    0 -0.75    0    0   10
[3,]    0 -0.75    0    0   30
[4,]    0 -0.75    0    1    5
[5,]    0 -0.75    0    1   10
[6,]    0 -0.75    0    1   30
 ```

#### runCell()
- applying the study procedure on a single cell (i.e., referred to as row in the `design` matrix above)


```r
single = runCell(shift_proportion = 50, 
                 shift_magnitude = .75, 
                 shift_type = 2,
                 parameters_type = 0,
                 test_length = 10)
```

#### runSelectedCells()
- applying the study procedure on a selected number of cell (i.e., it is a wrapper around `runCell()`)

```r
selected = runSelectedCells(design[c(1, 33, 99), ])
```

#### runSelectedCellsWithReplication()
- applying the study procedure on a selected number of cells and replicating the procedure `n` number of times (i.e., it is a wrapper around `runSelectedCells()`)

```r
selectedWithReplication = runSelectedCellsWithReplication(design[c(1, 33, 99), ], 100)
```

***Note:*** *Check the documentation for each function for more details regarding the inputs and the outputs.*

