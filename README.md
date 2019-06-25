<h1 align="center">reshiftsim</h1>

<p align="center">
    <a href="https://github.com/mihaiconstantin/reshiftsim/releases/latest">
        <img src="https://img.shields.io/badge/stable-v0.1.*-blue.svg?style=flat-square" alt="Latest Stable Release">
    </a>	
    <a href="https://opensource.org/licenses/MIT">
        <img src="https://img.shields.io/badge/license-MIT-yellow.svg?style=flat-square" alt="License">
    </a>
</p>

> Github repository name: **simulation-response-shift**

> Type: simulation package

## Authors

- [Mihai A. Constantin](https://constantinmihai.com) | Tilburg University
- *add collaborators here...*

## Description

**Data simulations** for a `Person-Fit Statistic` (l<sub>z</sub><sup>*</sup> for polychotomous items) using `R` and `C++`. Started as the First-Year Paper project at Tilburg University&mdash;Individual Differences and Assessment Research Master&mdash;under the supervison of dr. Wilco Emons.

Alongside the source files (`.cpp` and `.h`), this repository also contains the compiled files (`.dll` and `.o`) for a `64bit` Windows platform. You may also compile for your specific platform. Note that this package relies heavily on the [Rcpp](https://github.com/RcppCore/Rcpp) and [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo) packages. We advise you to have the appropriate tools ready (e.g., see [section 1.3](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf) from the `Rcpp` vignette) before attempting to build the package. We recommend using the [devtools](https://github.com/hadley/devtools) package for a straightforward installation and compilation.

## Important links

- this package is used in the manuscript within repository [`paper-response-shift`](https://github.com/mihaiconstantin/paper-response-shift), also available on OSF at https://osf.io/eg6yb

## Example

The code block below illustrates the main four functions in the package. Using them, one is able to replicate the simulations in the [manuscript](https://osf.io/eg6yb) or test how the procedure performs against a specific combination of factors (i.e., see the package documentation). Important information is also provided in the documentation associated with each function listed below.

### Exported functions

#### buildDesign()
- building a matrix based on the combinations of all the factors in the study
```r
# Specifing the factors.
shift.proportion = c(0, 30, 70)
shift.magnitude  = c(-.8, -.4, .4, .8)
shift.type       = c(0, 1, 2)
parameters.type  = c(0, 1)
test.length      = c(5, 15, 30)
	

# Building the factorial design.
design = buildDesign(shift.proportion,
		     shift.magnitude,
		     shift.type,
		     parameters.type,
		     test.length)
		     
> head(design)
     [,1]  [,2] [,3] [,4] [,5]
[1,]    0 -0.8    0    0    5
[2,]    0 -0.8    0    0   15
[3,]    0 -0.8    0    0   30
[4,]    0 -0.8    0    1    5
[5,]    0 -0.8    0    1   15
[6,]    0 -0.8    0    1   30
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

***Note:*** *Check the package documentation for each function for more details regarding the inputs and the outputs.*

## Installation

### Prerequisites

- see [section 1.3](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf) from the `Rcpp` vignette

### Installation steps

1. prepare all required build tools 
2. run `devtools::install_github("mihaiconstantin/simulation-response-shift")`

## Development

- feel free to email or send a PR request

## Roadmap

1. Features:
    - add new design factor to manipulate the response shift in the data matrix
2. Fixes:
    - &#x270B; check the code for errors and make sure that the `Rcpp` seeds work as expected
3. Improvements:
    - &#x2705; ~~update `README.md`~~

## License

The code in this repository is licensed under the [MIT license](https://opensource.org/licenses/MIT).
