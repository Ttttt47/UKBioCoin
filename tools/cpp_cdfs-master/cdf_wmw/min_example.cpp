#include <iostream>

//the cdf core code
#include "cdf_base.h"

//NOTE:
//may get a warning about undefining max/min; this can be fixed

int main(){

    //====================================//
    //EXAMPLE: cdf of normal distribution
    //====================================//
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "cdf of normal distribution" << std::endl;
    std::cout << "==========================" << std::endl;
    std::cout << std::endl;
    double x = 3.45;
    double mu = 1.2;
    double sigma = 1.5;
    //R code: pnorm(3.45, mean=1.2, sd=1.5)
    //Using R, the following values were obtained
    double r_p1 = 0.933192798731141914814202209527138620615005493164062500000000;
    double r_p1_log = -0.069143455612233992169102236857725074514746665954589843750000;
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


    //========================================//
    //EXAMPLE: cdf of Wilcoxon distribution 1
    //then do it with large sample, then do it with ties
    //========================================//
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "cdf of Wilcoxon distribution" << std::endl;
    std::cout << "============================" << std::endl;
    std::cout << std::endl;
    
    double q = 73;
    double m = 10;
    double n = 12;
    //R code: pnorm(3.45, mean=1.2, sd=1.5)
    //Using R, the following values were obtained
    r_p1 = 0.809514015396368247223790604039095342159271240234375000000000;
    r_p1_log = -0.211321192368411542306105843636032659560441970825195312500000;
    std::cout << "Computation of cdf of q=" << q << " for Wilcoxon distribution"; 
    std::cout << " with q equal to the sum of ranks for sample 1, " << std::endl;
    std::cout << " where sample 1 has m=" << m << " observations, and " << std::endl;
    std::cout << " sample 2 has n=" << n << "observations." << std::endl;
    std::cout << "Note: using the lower tail." << std::endl;

    std::cout << "(1) The p-value: " << r_p1 << std::endl;
    std::cout << "(2) The log of the p-value (computed directly): ";
    std::cout << r_p1_log << std::endl;

    std::cout << "(3) Check: Taking the log of (1): " << log(r_p1) << std::endl;
    std::cout << std::endl;

    prec = 50;
    std::cout << "Now computing with cdf_norm to " << prec;
    std::cout << " places." << std::endl;

    //now computing using cdf_base.h and cdf_base.cpp
    //the true is for "lower_tail"
    cpp_p1 = cdf_wilcoxon(q, m, n, true);
    cpp_p1_log = cdf_wilcoxon_log(q, m, n, true);
    
    //displaying results
    std::cout << "R:        "<< std::setprecision(prec) << r_p1 << std::endl;
    std::cout << "cdf_base: " << std::setprecision(prec) << cpp_p1 << std::endl;
    std::cout << std::endl;
    std::cout << "R (log):        "<< std::setprecision(prec) << r_p1_log;
    std::cout<< std::endl;
    std::cout << "cdf_base (log): " << std::setprecision(prec) << cpp_p1_log;
    std::cout << std::endl;




    return 0;
}
