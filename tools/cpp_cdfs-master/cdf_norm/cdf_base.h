#ifndef __CDF_NORM_H__
#define __CDF_NORM_H__

#include<float.h>
#include<math.h>
#include<string.h>
#include<ostream>
#include<iomanip>
#include "cdf_base.cpp"

//cdf of normal dist, no log
double cdf_norm(double x, double mu, double sigma, int lower_tail=true);

//cdf of normal dist, log
double cdf_norm_log(double x, double mu, double sigma, int lower_tail=true);

//quantile of normal dist, no log
double quantile_norm(double p, double mu, double sigma, int lower_tail=true);

//quantile of normal dist, log
double quantile_norm_log(double p, double mu, double sigma, int lower_tail=true);


#endif //__CDF_NORM_H__
