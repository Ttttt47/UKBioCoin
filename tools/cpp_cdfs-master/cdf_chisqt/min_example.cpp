#include <iostream>

//the cdf core code
#include "cdf_base.h"

//NOTE:
//may get a warning about undefining max/min; this can be fixed

int main(){

    //================================//
    //EXAMPLE: cdf of t-distribution
    //================================//
    double x = 0.7;
    double df = 2.0;
    //R code: pt(0.7, df=2.0, lower.tail=T, log.p=F)
    //Using R, the following values were obtained
    double r_p1 = 0.721803487683567279731278176768682897090911865234375000000000;
    double r_p1_log = -0.326002354859980858492463084985502064228057861328125000000000;
    std::cout << "Computation of cdf of x=" << x << " for t-distribution"; 
    std::cout << " with " << df << " degrees of freedom." << std::endl;
    std::cout << "Note: using the lower tail." << std::endl;

    std::cout << "(1) The p-value: " << r_p1 << std::endl;
    std::cout << "(2) The log of the p-value (computed directly): ";
    std::cout << r_p1_log << std::endl;

    std::cout << "(3) Check: Taking the log of (1): " << log(r_p1) << std::endl;
    std::cout << std::endl;

    int prec = 50;
    std::cout << "Now computing with cdf_t to " << prec;
    std::cout << " places." << std::endl;

    //now computing using cdf_base.h and cdf_base.cpp
    //the true is for "lower_tail"
    double cpp_p1 = cdf_t(x, df, true);
    double cpp_p1_log = cdf_t_log(x, df, true);
    
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
    //EXAMPLE: quantile of t-distribution
    //================================//
    std::setprecision(2);
    std::setprecision(2);
    double p = 0.95;
    df = 12.3;
    //R code: pt(0.7, df=2.0, lower.tail=T, log.p=F)
    //Using R, the following values were obtained
    double r_q1 = 1.778672546238235119275827855744864791631698608398437500000000;
    std::cout << "Computation of quantile of p=";
    std::cout << std::setprecision(3) << p << " for t-distribution"; 
    std::cout << " with " << df << " degrees of freedom." << std::endl;
    std::cout << "Note: using the lower tail." << std::endl;

    std::cout << "(1) The quantile-value: ";
    std::cout << std::setprecision(prec) << r_q1 << std::endl;

    prec = 50;
    std::cout << "Now computing with quantile_t to " << prec;
    std::cout << " places." << std::endl;

    //now computing using cdf_base.h and cdf_base.cpp
    //the true is for "lower_tail"
    double cpp_q1 = quantile_t(p, df);
    
    //displaying results
    std::cout << "R:        "<< std::setprecision(prec) << r_q1 << std::endl;
    std::cout << "cdf_base: " << std::setprecision(prec) << cpp_q1 << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    //======================================//
    //EXAMPLE cdf of chi-square distribution
    //======================================//

    //now computing for the chi-squared distribution
    std::setprecision(2);
    std::cout << "Computation of cdf of x=" << std::setprecision(2) <<  x;
    std::cout << " for chi-squared-"; 
    std::cout << "distribution with " << df << " degrees of freedom.";
    std::cout << std::endl;
    std::cout << "Computing with cdf_t to " << prec;
    std::cout << " places." << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    //R values
    x = 0.54;
    df = 2.0;
    //R code: pchisq(0.54, df=2.0, lower.tail=T, log.p=F)
    //Using R, the following values were obtained:
    double r_p2 = 0.236620505663146823982501132377365138381719589233398437500000;
    double r_p2_log = -1.441297663132672601804529222135897725820541381835937500000000;

    //cdf_base
    double cpp_p2 = cdf_chisq(x, df, true);
    double cpp_p2_log = cdf_chisq_log(x, df, true);


    //displaying results
    std::cout << "R:        "<< std::setprecision(prec) << r_p2 << std::endl;
    std::cout << "cdf_base: " << std::setprecision(prec) << cpp_p2 << std::endl;
    std::cout << std::endl;
    std::cout << "R (log):        "<< std::setprecision(prec) << r_p2_log;
    std::cout<< std::endl;
    std::cout << "cdf_base (log): " << std::setprecision(prec) << cpp_p2_log;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Note that cdf_base computes the log of the p-values directly";
    std::cout << " rather than first computing the 'raw' p-value and ";
    std::cout << "then taking the log." << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;


    //=============================================//
    //EXAMPLE: quantile of chi-square distribution
    //=============================================//
    p = 0.95;
    df = 12.3;
    //R code: pt(0.7, df=2.0, lower.tail=T, log.p=F)
    //Using R, the following values were obtained
    r_q1 = 21.428342879876034032804454909637570381164550781250000000000;
    std::cout << "Computation of quantile of p=";
    std::cout << std::setprecision(3) << p << " for chi-squared distribution"; 
    std::cout << " with " << df << " degrees of freedom." << std::endl;
    std::cout << "Note: using the lower tail." << std::endl;

    std::cout << "(1) The quantile-value: ";
    std::cout << std::setprecision(prec) << r_q1 << std::endl;

    prec = 50;
    std::cout << "Now computing with quantile_t to " << prec;
    std::cout << " places." << std::endl;

    //now computing using cdf_base.h and cdf_base.cpp
    //the true is for "lower_tail"
    cpp_q1 = quantile_chisq(p, df);
    
    //displaying results
    std::cout << "R:        "<< std::setprecision(prec) << r_q1 << std::endl;
    std::cout << "cdf_base: " << std::setprecision(prec) << cpp_q1 << std::endl;
    std::cout << std::endl;


    return 0;
}
