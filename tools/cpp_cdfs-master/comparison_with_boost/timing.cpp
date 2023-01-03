#include <iostream>
#include <cstdio>
#include <ctime>

// the boost libraries
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
//the cdf core code
#include "cdf_base.h"

int main() {
    std::clock_t start;
    double boost_duration_t;
    double cdf_base_duration_t;
    double cdf_base_duration_t_log;

    double boost_duration_chisq;
    double cdf_base_duration_chisq;
    double cdf_base_duration_chisq_log;

    //a few random choices
    double x1 = 0.7;
    double df1 = 2.0;

    double x2 = 1.2;
    double df2 = 3.5;

    double x3 = 6.4;
    double df3 = 7.55;

    double x4 = 12.5;
    double df4 = 20.06;

    double x5 = 0.04;
    double df5 = 1.23;


    //----------------------------------------------------------------//
    //start clock for boost t-dist
    start = std::clock();
    
    //timing for boost
    boost::math::students_t t_dist1(df1);
    double boost_p1 = cdf(t_dist1, fabs(x1));

    boost::math::students_t t_dist2(df2);
    double boost_p2 = cdf(t_dist2, fabs(x2));

    boost::math::students_t t_dist3(df3);
    double boost_p3 = cdf(t_dist3, fabs(x3));

    boost::math::students_t t_dist4(df4);
    double boost_p4 = cdf(t_dist4, fabs(x4));

    boost::math::students_t t_dist5(df5);
    double boost_p5 = cdf(t_dist5, fabs(x5));

    //end clock
    boost_duration_t = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //----------------------------------------------------------------//


    //----------------------------------------------------------------//
    //timing for cdf_base t
    //start clock
    start = std::clock();

    double cpp_p1 = cdf_t(x1, df1, true);
    double cpp_p2 = cdf_t(x2, df2, true);
    double cpp_p3 = cdf_t(x3, df3, true);
    double cpp_p4 = cdf_t(x4, df4, true);
    double cpp_p5 = cdf_t(x5, df5, true);

    //end clock
    cdf_base_duration_t = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //----------------------------------------------------------------//


    //----------------------------------------------------------------//
    //timing for cdf_base t log
    //start clock
    start = std::clock();

    double cpp_p1_log = cdf_t_log(x1, df1, true);
    double cpp_p2_log = cdf_t_log(x2, df2, true);
    double cpp_p3_log = cdf_t_log(x3, df3, true);
    double cpp_p4_log = cdf_t_log(x4, df4, true);
    double cpp_p5_log = cdf_t_log(x5, df5, true);

    //end clock
    cdf_base_duration_t_log = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //----------------------------------------------------------------//
   

    //----------------------------------------------------------------//
    //timing for boost_chisq
    //start clock
    start = std::clock();

    boost::math::chi_squared chisq_dist1(df1);
    double boost_chi1 = cdf(chisq_dist1, x1);

    boost::math::chi_squared chisq_dist2(df2);
    double boost_chi2 = cdf(chisq_dist2, x2);

    boost::math::chi_squared chisq_dist3(df3);
    double boost_chi3 = cdf(chisq_dist3, x3);

    boost::math::chi_squared chisq_dist4(df4);
    double boost_chi4 = cdf(chisq_dist4, x4);

    boost::math::chi_squared chisq_dist5(df5);
    double boost_chi5 = cdf(chisq_dist5, x5);

    //end clock
    boost_duration_chisq = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //----------------------------------------------------------------//
    

    //----------------------------------------------------------------//
    //timing for cdf_base_chisq 
    //start clock
    start = std::clock();

    double cpp_chi_p1 = cdf_chisq(x1, df1, true);
    double cpp_chi_p2 = cdf_chisq(x2, df2, true);
    double cpp_chi_p3 = cdf_chisq(x3, df3, true);
    double cpp_chi_p4 = cdf_chisq(x4, df4, true);
    double cpp_chi_p5 = cdf_chisq(x5, df5, true);

    //end clock
    cdf_base_duration_chisq = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //----------------------------------------------------------------//


    //----------------------------------------------------------------//
    //timing for cdf_base_chisq 
    //start clock
    start = std::clock();

    double cpp_chi_p1_log = cdf_chisq_log(x1, df1, true);
    double cpp_chi_p2_log = cdf_chisq_log(x2, df2, true);
    double cpp_chi_p3_log = cdf_chisq_log(x3, df3, true);
    double cpp_chi_p4_log = cdf_chisq_log(x4, df4, true);
    double cpp_chi_p5_log = cdf_chisq_log(x5, df5, true);

    //end clock
    cdf_base_duration_chisq_log = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //----------------------------------------------------------------//
    std::cout << "t: boost duration: "<< boost_duration_t << std::endl;
    std::cout << "t: cdf_base duration: "<< cdf_base_duration_t << std::endl;
    std::cout << "t: cdf_base duration log: "<< cdf_base_duration_t_log << std::endl;

    std::cout << "chisq: boost duration: "<< boost_duration_chisq << std::endl;
    std::cout << "chisq: cdf_base duration: "<< cdf_base_duration_chisq << std::endl;
    std::cout << "chisq: cdf_base duration log: "<< cdf_base_duration_chisq_log << std::endl;


}
