
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tvRRR –

## Time-varying reduced rank regression using state-space models

<!-- badges: start -->
<!-- badges: end -->

The tvRRR package implements time-varying reduced rank regression as
proposed in *A state-space approach to time-varying reduced rank
regression* (Brune, Bura and Scherrer, 2021+).

## Installation

The development version can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("b-brune/tvRRR")
```

## Model

The model that is fitted is a reduced-rank regression with time-varying
coefficient matrices. We implement two types of time-variation, called
model (A) and model (B):

(*A*)  *y*<sub>*t*</sub> = *α*<sub>*t*</sub>*β*′*x*<sub>*t*</sub> + *Γ**u*<sub>*t*</sub> + *ε*<sub>*t*</sub>

(*B*)  *y*<sub>*t*</sub> = *α**β*<sub>*t*</sub>′*x*<sub>*t*</sub> + *Γ**u*<sub>*t*</sub> + *ε*<sub>*t*</sub>

Model fitting is carried out using a Gaussian linear state-space model
and the Kalman filter.

## Example application

The function `dataset()` draws data from models (A) and (B) with
different specifications of the time-variation in the parameter matrix.

``` r
library(tvRRR)

set.seed(712)
dat <- dataset("VARbreak", p = 5, d = 2, q = 5, t = 100, model = "A")
X <- dat$X; y <- dat$y
```

The main function of the dataset is `tvRRR()`. It automatically selects
the rank of the model using BIC (in case it is unknown) and fits the
tvRRR model.

``` r
fit <- tvRRR(X = X, y = y, select_rank = TRUE, d_max = 5, silent = TRUE, 
             model = "A")

fit
#> Time-varying reduced rank regression model of type A with latent rank d = 2 
#> 
#> Convergence information: 
#> 
#> Model selected based on BIC
#> 
#> Tried d = 1 to d = 5
#> Rank selected: d = 2 with BIC 558.277561
#> 
#> Model diagnostics / Convergence information 
#> 
#> Algorithm stopped after 29 iterations. 
#>  Relative likelihood difference: 0.000378
#>  Difference in Omega: 0.000968
#>  Difference in beta: 0.00052
#>  Difference in Sigma: 0.000436 
#> 
#> Estimated state covariance matrix: 
#> 
#>              [,1]         [,2]
#> [1,]  0.047435252 -0.003688764
#> [2,] -0.003688764  0.086516325
#> 
#>  Proportion of variance explained per time series: 
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.4795375 0.5628192 0.7748703 0.5245667 0.4412078
```

`tvRRR()` returns a list of the filtered time-varying parameters
(one-step-ahead, filtered and smoothed), as well as the fitted
parameters and some information on the data log-likelihood and
convergence.

For more information on the model fitting algorithm see the paper, the
accompanying vignette, and the function’s documentation.
