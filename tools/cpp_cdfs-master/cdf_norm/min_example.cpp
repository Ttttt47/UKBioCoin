#include <iostream>

#include "cdf_base.h"
//may get a warning about undefining max/min; this can be fixed

int main(){

    //====================================//
    //EXAMPLE: cdf of normal distribution
    //====================================//
    
    double x = 3.45;
    double mu = 1.2;
    double sigma = 1.5;
    //R code: pnorm(3.45, mean=1.2, sd=1.5)
    //Using R, the following values were obtained
    double r_p1 = 0.933192798731141914814202209527138620615005493164062500000000;
    double r_p1_log = -0.069143455612233992169102236857725074514746665954589843750000;
    std::cout << std::endl;
    std::cout << "Computation of cdf of x=" << x << " for normal distribution"; 
    std::cout << " with mean " << mu << " and sd " << sigma << "." << std::endl;
    std::cout << "Note: using the lower tail." << std::endl;

    std::cout << "(1) The p-value: " << r_p1 << std::endl;
    std::cout << "(2) The log of the p-value (computed directly): ";
    std::cout << r_p1_log << std::endl;

    std::cout << "(3) Check: Taking the log of (1): " << log(r_p1) << std::endl;
    std::cout << std::endl;

    int prec = 50;
    std::cout << "Now computing with cdf_norm to " << prec;
    std::cout << " places." << std::endl;

    //now computing using cdf_base.h and cdf_base.cpp
    //the true is for "lower_tail"
    double cpp_p1 = cdf_norm(x, mu, sigma, true);
    double cpp_p1_log = cdf_norm_log(x, mu, sigma, true);
    
    //displaying results
    std::cout << "R:        "<< std::setprecision(prec) << r_p1 << std::endl;
    std::cout << "cdf_base: " << std::setprecision(prec) << cpp_p1 << std::endl;
    std::cout << std::endl;
    std::cout << "R (log):        "<< std::setprecision(prec) << r_p1_log;
    std::cout<< std::endl;
    std::cout << "cdf_base (log): " << std::setprecision(prec) << cpp_p1_log;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;


    //================================//
    //EXAMPLE: quantile of normal distribution
    //================================//
    std::setprecision(2);
    std::setprecision(2);
    double p = 0.975;
    mu = 0;
    sigma = 1;
    //R code:
    //> q <- qnorm(0.975, 0, 1)
    //> sprintf("%.60f", q)
    //[1] "1.959963984540053605343246090342290699481964111328125000000000"
    //Using R, the following values were obtained
    double r_q1 = 1.959963984540053605343246090342290699481964111328125000000000;
    std::cout << "Computation of quantile of p=";
    std::cout << std::setprecision(3) << p << " for t-distribution"; 
    std::cout << " with mean " << mu << " and sd " << sigma << "." << std::endl;
    std::cout << "Note: using the lower tail." << std::endl;

    std::cout << "(1) The quantile-value: ";
    std::cout << std::setprecision(prec) << r_q1 << std::endl;

    prec = 50;
    std::cout << "Now computing with quantile_t to " << prec;
    std::cout << " places." << std::endl;

    //now computing using cdf_base.h and cdf_base.cpp
    //the true is for "lower_tail"
    double cpp_q1 = quantile_norm(p, mu, sigma);
    
    //displaying results
    std::cout << "R:        "<< std::setprecision(prec) << r_q1 << std::endl;
    std::cout << "cdf_base: " << std::setprecision(prec) << cpp_q1 << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    return 0;
}
