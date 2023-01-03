# Comparison with Boost

**Note: this was the first version of the repo, where I was interested in 
comparing the performance of R's implementations with Boost's.**

This code is taken from the R core code for computing the cdfs of the t- and chiq-squared distributions. All of the code was created by the R developers, was originally written in C, and there are minor modifications to port it to C++. Please note that there is no warranty for the code! A minimal example suggests it gives the same answers as the Boost library.


 * Using the R code code.
 * Compared to Boost implementation

The key files are:
 * `cdf_base.cpp`
 * `cdf_base.h`

Include those two files in your directory, and the line
```
#include "cdf_base.h"
```

Please see `min_example.cpp` for details.


At the moment there is only a basic runtime check, but it seems computing the log of the p-value for the t-distribution cdf is an order of magnitude faster than the Boost computation of the 'raw' p-value.

 ```  
 t: boost duration: 6.3e-05
 t: cdf_base duration: 4.6e-05
 t: cdf_base duration log: 7e-06
 chisq: boost duration: 9e-06
 chisq: cdf_base duration: 1.9e-05
 chisq: cdf_base duration log: 6e-06
 ```  

Please see `min_example.cpp` or `minmin_example.cpp` for usage of the `cdf_base.h` functions.

Make the `.sh` files executable in order to run the example, tests, etc.
```
chmod a+x run_minmin_example.sh
./run_minmin_example.sh
```

Note that the `doctest.h` file used in the tests is from [https://github.com/onqtam/doctest](https://github.com/onqtam/doctest)

