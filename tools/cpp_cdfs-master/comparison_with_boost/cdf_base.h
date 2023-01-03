#include "cdf_base.cpp"

//Cdf which computes the p-value of the t distribution
//
//double cdf_t(double x, double n, int lower_tail){
//x: the quantile value
//n: the degrees of freedom
//lower_tail: a boolean of whether or not to use the lower tail
double cdf_t(double, double, int);

//Cdf which computes the log of the p-value
//cdf_t_log: same parameters as for cdf_t
double cdf_t_log(double, double, int);

//Cdf which computes the p-value of the chi-squared distribution
//x: the quantile value
//df: the degrees of freedom
//lower_tail: a boolean of whether or not to use the lower tail
double cdf_chisq(double, double, int);

//Cdf which computes the log of the p-value
//cdf_chisq_log: same parameters as for cdf_chisq
double cdf_chisq_log(double, double, int);
