#ifndef __CDFS_H__
#define __CDFS_H__

#include<float.h>
#include<math.h>
#include<string.h>
#include<ostream>
#include<iomanip>
#include <iostream>
#include "cdf_base.cpp"


//Cdf which computes the p-value of the t distribution
//
//
//double cdf_t(double x, double n, int lower_tail){
//x: the quantile value
//n: the degrees of freedom
//lower_tail: a boolean of whether or not to use the lower tail
double cdf_t(double q, double df, int lower_tail=true);

//Cdf which computes the log of the p-value
//cdf_t_log: same parameters as for cdf_t
double cdf_t_log(double q, double df, int lower_tail=true);

//Cdf which computes the p-value of the chi-squared distribution
//x: the quantile value
//df: the degrees of freedom
//lower_tail: a boolean of whether or not to use the lower tail
double cdf_chisq(double q, double df, int lower_tail=true);

//Cdf which computes the log of the p-value
//cdf_chisq_log: same parameters as for cdf_chisq
double cdf_chisq_log(double q, double df, int lower_tail=true);

//Quantile function for Student's t-distribution
double quantile_t(double p, double df, int lower_tail=true);

//Quantile function for chi-squared distribution
double quantile_chisq(double p, double df, int lower_tail=true);


//cdf of normal dist, no log
double cdf_normal(double x, double mu, double sigma, int lower_tail=true);

//cdf of normal dist, log
double cdf_normal_log(double x, double mu, double sigma, int lower_tail=true);

//quantile of normal dist, no log
double quantile_normal(double p, double mu, double sigma, int lower_tail=true);

//quantile of normal dist, log
double quantile_normal_log(double p, double mu, double sigma, int lower_tail=true);


//cdf of wilcoxon, no log
double cdf_wilcoxon(double q, double m, double n, int lower_tail=true);


//cdf of wilcoxon, no log
double cdf_wilcoxon_log(double q, double m, double n, int lower_tail=true);


//cdf of Gamma distribution, log_p=false
double cdf_gamma(double x, double shape, double scale, int lower_tail=true);


//cdf of Gamma distribution, log_p=true
double cdf_gamma_log(double x, double shape, double scale, int lower_tail=true);


//quantile function of Gamma distribution, log_p=false
double quantile_gamma(double p, double shape, double scale, int lower_tail=true);

//quantile function of Gamma distribution, log_p=true
double quantile_gamma_log(double p, double shape, double scale, int lower_tail=true);


//noncentral chisq cdf, no log 
double cdf_noncentral_chisq(double x, double df, double ncp=0, int lower_tail=true);


//noncentral chisq cdf, log
double cdf_noncentral_chisq_log(double x, double df, double ncp=0, int lower_tail=true);


//noncentral chi-squared distribution quantile function
double quantile_noncentral_chisq(double p, double df, double ncp=0, int lower_tail=true);


//beta cdf, no log 
double cdf_beta(double x, double a, double b, int lower_tail=true);

//beta cdf, with log 
double cdf_beta_log(double x, double a, double b, int lower_tail=true);


// poisson cdf, no log
double cdf_poisson(double x, double lambda, int lower_tail=true);
   

// poisson cdf, with log
double cdf_poisson_log(double x, double lambda, int lower_tail=true);


// binomial cdf, no log
double cdf_binomial(double x, double n, double p, int lower_tail=true);


// binomial cdf, with log
double cdf_binomial_log(double x, double n, double p, int lower_tail=true);

// F distribution cdf, no log
double cdf_f(double x, double df1, double df2, int lower_tail=true);

// F distribution cdf, no log
double cdf_f_log(double x, double df1, double df2, int lower_tail=true);


// Exponential distribution cdf, no log
double cdf_exp(double x, double rate, int lower_tail=true);


// Exponential distribution cdf, with log
double cdf_exp_log(double x, double rate, int lower_tail=true);


// Geometric distribution cdf, no log
double cdf_geom(double x, double p, int lower_tail=true);


// Geometric distribution cdf, with log
double cdf_geom_log(double x, double p, int lower_tail=true);


// Cauchy distribution cdf, no log
double cdf_cauchy(double x, double loc, double scale, int lower_tail=true);


// Cauchy distribution cdf, with log
double cdf_cauchy_log(double x, double loc, double scale, int lower_tail=true);


// Weibull distribution cdf, no log
double cdf_weibull(double x, double shape, double scale, int lower_tail);


// Weibull distribution cdf, with log
double cdf_weibull_log(double x, double shape, double scale, int lower_tail);


//hypergeometric distribution cdf, no log
//                   q        NR          NB         ND   
/**                  q         m           n         k   
 *               quantile   #red      #black        #drawn */
double cdf_hypergeometric(double x, double m, double n, double k, int lower_tail=true);


double cdf_hypergeometric_log (double x, double m, double n, double k, int lower_tail=true);


//cdf of lognormal dist, no log
double cdf_lognormal(double x, double meanlog, double sdlog, int lower_tail);

//cdf of lognormal dist, log
double cdf_lognormal_log(double x, double meanlog, double sdlog, int lower_tail);

#endif //__CDFS_H__
