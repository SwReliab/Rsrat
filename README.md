
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rsrat

Rsrat provides the package to evalute the software reliability from the
fault data collected in the testing phase. Rsrat can use two types of
data; fault-detection time data and its grouped data. The
fault-detection time data is a sequence of time intervals of fault
detection times (CPU time, etc). Also its grouped data is a sequence of
the number of detected faults for each time interval (per a working day,
per a week, etc). The reliability evaluation is based on the software
reliability growth model with NHPP (non-homogeneous Poisson process).

## Installation

You can install Rsrat from github with:

``` r
# install.packages("devtools")
devtools::install_github("okamumu/Rsrat")
```

## Example

This is an example of the estimation of software reliability growth
models from a fault data (tohma).

``` r
### load library
library(Rsrat)

### load example data
data(dacs)

### tohma is a grouped data
tohma
#>   [1]  5  5  5  5  6  8  2  7  4  2 31  4 24 49 14 12  8  9  4  7  6  9  4
#>  [24]  4  2  4  3  9  2  5  4  1  4  3  6 13 19 15  7 15 21  8  6 20 10  3
#>  [47]  3  8  5  1  2  2  2  7  2  0  2  3  2  7  3  0  1  0  1  0  0  1  1
#>  [70]  0  0  1  1  0  0  0  1  2  0  1  0  0  0  0  0  0  2  0  0  0  0  0
#>  [93]  0  0  0  1  0  0  0  1  0  0  1  0  0  1  0  0  1  0  1

### Esimate all models and select the best one in terms of AIC
(result <- fit.srm.nhpp(fault=tohma))
#> Model name: LXVMinSRM
#>    omega    loclog  scalelog  
#> 481.7029   -3.4642    0.6637  
#> Maximum LLF: -316.2599 
#> AIC: 638.5198 
#> Convergence: TRUE

### Draw the graph 
mvfplot(fault=tohma, mvf=list(result$srm))
```

<img src="man/figures/README-example1-1.png" width="100%" />

The second example illustrates the estimation for two specified models.

``` r
### All models in the package
srm.models
#>  [1] "exp"    "gamma"  "pareto" "tnorm"  "lnorm"  "tlogis" "llogis"
#>  [8] "txvmax" "lxvmax" "txvmin" "lxvmin"

### Estimate two models and no select
(result <- fit.srm.nhpp(fault=tohma, srm.names=c("exp", "gamma"), selection=NULL))
#> $exp
#> Model name: ExpSRM
#>    omega      rate  
#> 497.2912    0.0308  
#> Maximum LLF: -359.8777 
#> AIC: 723.7555 
#> Convergence: TRUE 
#> 
#> 
#> $gamma
#> Model name: GammaSRM
#>     omega      shape       rate  
#> 483.52301    1.88475    0.06447  
#> Maximum LLF: -319.5695 
#> AIC: 645.139 
#> Convergence: TRUE

### Draw the graph
mvfplot(fault=tohma, mvf=lapply(result, function(m) m$srm))
```

<img src="man/figures/README-example2-1.png" width="100%" />

The third example shows the case where the fault data are fault
detection data.

``` r
### fault-detection time data
#### Time intervals for all faults
#### The last value is a negative value, that indicates the time interval in which there is no fault detection after the last fault detection.
sys1
#>   [1]     3    30   113    81   115     9     2    91   112    15   138
#>  [12]    50    77    24   108    88   670   120    26   114   325    55
#>  [23]   242    68   422   180    10  1146   600    15    36     4     0
#>  [34]     8   227    65   176    58   457   300    97   263   452   255
#>  [45]   197   193     6    79   816  1351   148    21   233   134   357
#>  [56]   193   236    31   369   748     0   232   330   365  1222   543
#>  [67]    10    16   529   379    44   129   810   290   300   529   281
#>  [78]   160   828  1011   445   296  1755  1064  1783   860   983   707
#>  [89]    33   868   724  2323  2930  1461   843    12   261  1800   865
#> [100]  1435    30   143   108     0  3110  1247   943   700   875   245
#> [111]   729  1897   447   386   446   122   990   948  1082    22    75
#> [122]   482  5509   100    10  1071   371   790  6150  3321  1045   648
#> [133]  5485  1160  1864  4116 -2526

### Esimate
(result <- fit.srm.nhpp(time=sys1[sys1>=0], te=-sys1[sys1<0]))
#> Warning in emfit(srm, data, initialize = TRUE, maxiter = con$maxiter,
#> reltol = con$reltol, : Did not converge to MLE by max iteration.
#> Model name: LXVMinSRM
#>    omega    loclog  scalelog  
#>  165.955   -10.640     1.453  
#> Maximum LLF: 299.439 
#> AIC: -592.8779 
#> Convergence: TRUE

### Draw the graph
mvfplot(time=sys1[sys1>=0], te=-sys1[sys1<0], mvf=list(result$srm))
```

<img src="man/figures/README-example3-1.png" width="100%" />
