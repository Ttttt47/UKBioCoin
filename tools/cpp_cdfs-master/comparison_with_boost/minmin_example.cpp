#include <iostream>

// the boost libraries
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
//the cdf core code
#include "cdf_base.h"

//NOTE:
//may get a warning about undefining max/min; this can be fixed

int main(){
    //EXAMPLE for t-distribution
    //
    double x = 0.7;
    double df = 2.0;
    //R code: pt(0.7, df=2.0, lower.tail=T, log.p=F)
    //Using R, the following values were obtained
    double r_p1 = 0.721803487683567279731278176768682897090911865234375000000000;
    double r_p1_log = -0.326002354859980858492463084985502064228057861328125000000000;

    //the true is for "lower_tail"
    double cpp_p1 = cdf_t(x, df, true);
    double cpp_p1_log = cdf_t_log(x, df, true);
    
    //now computing using boost
    boost::math::students_t t_dist(df);
    double boost_p1 = cdf(t_dist, fabs(x));

    std::cout << cpp_p1 << std::endl;
    std::cout << boost_p1 << std::endl;
    std::cout << cpp_p1_log << std::endl;
    std::cout << log(boost_p1) << std::endl;
}
