# C++ implementations of cdfs from R source

Implementation of the following cumulative distribution functions from
the R source.

  - Normal distribution (cdf and quantile functon)
  - Student t-distribution (cdf and quantile function)
  - Gamma distribution (cdf and quantile function)
  - Chi-square distribution (cdf and quantile function)
  - Non-central chi-square distribution (cdf and quantile function)
  - Beta distribution
  - Poisson distribution cdf
  - Binomial distribution cdf
  - F distribution cdf
  - Wilcoxon distribution cdf 
  - Exponential
  - Geometric
  - Cauchy
  - Weibull
  - Hypergeometric
  - Lognormal


## Example 


Here is an example of the gamma distribution cdf being called (not the shape/scale 
parametrisation).

```
#include <iostream>

//the cdf core code, make sure the path is correct!
#include "cdf_base.h"

int main(){
    double x = 0.12;
    double shape = 1.2;
    double scale = 2.3;
    //Using R, the following values were obtained:
    double r_p = 0.0254992146698800326021228102035820484161376953125000000000;

    //cdf_base
    double p = cdf_gamma(x, shape, scale, true);

    //displaying results
    int prec=12;
    std::cout << "cdf_base: " << std::setprecision(prec) << p << std::endl;
    std::cout << "R:        " << std::setprecision(prec) << r_p << std::endl;
    return 0;
}
```

Please see `min_example.cpp` for usage of the `cdf_base.h` functions, or look
as `tests/test_h.cpp`, to see how to call other functions; generally follows the 
R syntax.
