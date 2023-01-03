# cpp\_cdfs: Cumulative Distribution Functions in C++

This repository contains code for cumulative distribution functions (cdfs)
of various probability distributions. **Almost all of the code is from the
source code of the R programming language and was created by the R Core Team.**
That code (and therefore this repository) is freely available under a
[GNU General Public License](https://www.R-project.org/Licenses/).

Here is a link to the [R programming language](https://www.r-project.org)
and a link to the GitHub repository for
the [R source code](https://github.com/wch/r-source) from which this
repo was created.

If you want to read about my motivation for creating and sharing this repo,
see below.



## A guide to the subrepositories

As shown in the table below, The various cdfs are distributed across several
subrepos, each of which contains a header file `cdf_base.h` exposing the
relevant functions. Each subrepo has its own page with example calls.

One of the repos contains all the cdfs (and is significantly longer, with
its `cdf_base.cpp` over 6000 lines, while the corresponding file for the
normal distribution is under 1000 lines). Most of the 'effort' in the full
file goes towards computing the beta and gamma cdfs.

|  Subfolder  |           |Distributions|        |             |           |
| :---------  | :-------: | :---------: | :----: | :---------: | :--------:|
|             | Normal    | Student's t | Gamma  | Chi-squared | Wilcoxon  |
| cdf_norm    | &check;   |             |        |             |           |
| cdf_chisqt  |           | &check;     | &check;| &check;     |           |
| cdf_wmw     | &check;   |             |        |             | &check;   |
| cdfs        | &check;   | &check;     | &check;| &check;     | &check;   |

In the full `cdfs` subfolder, there are also cdfs for the following distributions:

  - Non-central chi-squared
  - Beta
  - Poisson
  - Binomial
  - F distribution
  - Exponential
  - Geometric
  - Cauchy
  - Weibull
  - Hypergeometric
  - Lognormal distribution

Note that the Gamma cdf uses the shape/scale parametrisation.

## So why make this repository?

When writing code for statistical applications, the final result may end up
being a test statistic which needs to be 'turned into' a *p*-value.

In one of my research papers I needed the cdf for the chi-squared
distribution. I thought it would be a simple matter
to get the source code from `pchisq.c` in the R source code, and then call that
function. However, the `pchisq` function called the `pgamma` function
(unsurprisingly), so the code from `pgamma.c` was also needed. But then
one of the subroutines required `dnorm`, and then there were references to
various constants not defined in those files...in other words, there was a large
amount of interdependency between the various routines.

So, the idea was to create a single `.cpp` file that had all of the variables
and functions necessary to call the relevant cdf. This involved a lot
of `grep`ing to track down the various
variables/constants/functions; see `checklist.txt` if interested.
(Not particularly hard, just took a while; hope it save others some time.)


I then saw
a [Stackoverflow](https://stackoverflow.com/questions/55010268/is-there-a-built-in-chi-square-cdf-function-in-c) post about someone looking for the implementation
for the cdf of the noncentral chi-square distribution, and decided to share
this.

### Why not just use R?

R is probably the best language for quickly performing a statistical analysis
on a data set, but there are parts of it that are a bit slow (e.g. for loops).

One idea is to write the code for the computationally-intensive parts in C++,
and then call these bits from R using `Rcpp`, which is
a [great package](https://cran.r-project.org/web/packages/Rcpp/index.html).
One can also call R functions (e.g. the cdfs) from
C++ using `Rcpp`. I use `Rcpp` and highly recommend it.

Lately, however, I have decided to make my statistical software packages
available in both R and Python, in order to serve both communities. So I now
write the core code in C++ and then wrap it using `Rcpp` for R, and use
Cython to wrap it
for Python. However, at some point in my methods I usually need to call a cdf
function, which is why I need the cdfs in 'pure C++'. I then started to seeing
if I could use parts of the R source code for the cdfs, which led to this repo.

I have **huge respect** for the team of programmers that wrote the underlying
C code and designed the subroutines; to me, it is absolutely amazing.


## Notes

### Minimal working examples

Minimal working examples are provided in each subrepo; using GCC is installed,
from a command line simply run `./run_minexample.sh`.

Or, check the code in `minexample.cpp` in order to see the function calls.

There are a few cdfs that are not called in the `minexample.cpp` file - in that
case, check the tests in the `tests` subfolder; each cdf is tested at least
once.


### Unit tests

There are unit tests in the `tests` subfolder which calls runs
a few tests using the awesome and lightweight
[doctest](https://github.com/doctest/doctest) testing suite. Note that only one
copy of the `doctest.h` file is stored in the `doctest` subfolder; so if you
want to run the rests yourself, you will need that file (and may need to change
the path in `tests/test_h.cpp` if you do not download the whole repo).


<!--
### Dependencies

If you are interested in seeing how the implementations of the distributions depend
on each other, see here.
-->


### Error handling

Occasionally the core C code in the R source would run error checks, e.g.
`R_CheckUserInterrupt();`.

I tried hard to incorporate this, but ended up getting stuck; I could not
get the `SEXP` type to be defined and kept getting compile errors so eventually
gave up. In place of these checks, I just made a call to the `warning` function
which would print a warning to screen. It is not ideal, if anyone can fix this
in a better way, please get in touch/submit a pull request.
