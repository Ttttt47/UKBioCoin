#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../../doctest/doctest.h"
// or if you include doctest.h in the same folder:
// #include "../../doctest/doctest.h"

//include others
#include <iostream>
#include "../cdf_base.h"


TEST_CASE("cdf norm values"){
    double x = 3.45;
    double mu = 1.2;
    double sigma = 1.5;

    // p <- pnorm(3.45, mean=1.2, sd=1.5)
    // sprintf("%.60f", p)
    // "0.933192798731141914814202209527138620615005493164062500000000"
    //
    // p <- pnorm(3.45, mean=1.2, sd=1.5, lower.tail=F, log.p=F)
    // sprintf("%.60f", p)
    // "0.066807201268858071308009982658404624089598655700683593750000"
    //
    // p <- pnorm(3.45, mean=1.2, sd=1.5, lower.tail=T, log.p=T)
    // sprintf("%.60f", p)
    // "-0.069143455612233992169102236857725074514746665954589843750000"
    //
    // p <- pnorm(3.45, mean=1.2, sd=1.5, lower.tail=F, log.p=T)
    // sprintf("%.60f", p)
    // "-2.705944400823889761653617824777029454708099365234375000000000"


    double r_p1 = 0.933192798731141914814202209527138620615005493164062500000000;
    double r_p1_up = 0.066807201268858071308009982658404624089598655700683593750000;
    double r_p1_log = -0.069143455612233992169102236857725074514746665954589843750000;
    double r_p1_up_log= -2.705944400823889761653617824777029454708099365234375000000000;

    double e_cpp_p1 = cdf_norm(x, mu, sigma);
    double e_cpp_p1_log = cdf_norm_log(x, mu, sigma,  true);
    double e_cpp_p1_up = cdf_norm(x, mu, sigma, false);
    double e_cpp_p1_up_log = cdf_norm_log(x, mu, sigma, false);

    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}



TEST_CASE("quantile norm values"){
    double p = 0.95;
    double mu = 0;
    double sigma = 1;
    // R code
    // > q <- qnorm(0.95, 0, 1)
    // > sprintf("%.60f", q)
    // [1] "1.644853626951471525785564153920859098434448242187500000000000"
    double soln = 1.644853626951471525785564153920859098434448242187500000000000;
    double ans = quantile_norm(p=p, mu=mu, sigma=sigma);

    double eps = 1e-50;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));


     p = 0.01;
     mu = 0;
     sigma = 1;
     // > q <- qnorm(0.01, 0, 1)
     // > sprintf("%.60f", q)
     // [1] "-2.326347874040840757459136511897668242454528808593750000000000"
     soln = -2.326347874040840757459136511897668242454528808593750000000000;
     ans = quantile_norm(p, mu, sigma);

     eps = 1e-50;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));

}
