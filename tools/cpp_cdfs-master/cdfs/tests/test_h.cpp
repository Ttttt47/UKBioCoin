#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../../doctest/doctest.h"

//include others
#include <iostream>
#include "../cdf_base.h"

//=========================================================================//

TEST_CASE("cdf t values"){

    //x <- pt(0.7, df=2.0, lower.tail=T, log.p=F)
    // sprintf("%.60f", x)
    // [1] "0.721803487683567279731278176768682897090911865234375000000000"

    // x <- pt(0.7, df=2.0, lower.tail=T, log.p=T)
    // sprintf("%.60f", x)
    // [1] "-0.326002354859980858492463084985502064228057861328125000000000"

    // x <- pt(0.7, df=2.0, lower.tail=F, log.p=F)
    // sprintf("%.60f", x)
    // [1] "0.278196512316432720268721823231317102909088134765625000000000"

    // x <- pt(0.7, df=2.0, lower.tail=F, log.p=T)
    // sprintf("%.60f", x)
    // [1] "-1.279427536176580559867943520657718181610107421875000000000000"


    double x = 0.7;
    double df = 2.0;
    double r_p1 = 0.721803487683567279731278176768682897090911865234375000000000;
    double r_p1_log = -0.326002354859980858492463084985502064228057861328125000000000;
    double r_p1_up = 0.278196512316432720268721823231317102909088134765625000000000;
    double r_p1_up_log= -1.279427536176580559867943520657718181610107421875000000000000;


    double e_cpp_p1 = cdf_t(x, df, true);
    double e_cpp_p1_log = cdf_t_log(x, df, true);
    double e_cpp_p1_up = cdf_t(x, df, false);
    double e_cpp_p1_up_log = cdf_t_log(x, df, false);
    
    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}//testing pt



TEST_CASE("quantile t values"){
    double p = 0.95;
    double df = 12.3;
    // q <- qt(0.95, df=12.3)
    // sprintf("%.60f", q)
    // "1.778672546238235119275827855744864791631698608398437500000000"
    double soln = 1.778672546238235119275827855744864791631698608398437500000000;
    double ans = quantile_t(p=p, df=df);

    double eps = 1e-50;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));


     p = 0.01;
     df = 3.45;
     // q <- qt(0.01, 3.45)
     // sprintf("%.60f", q)
     // "-4.099527305346944316966073529329150915145874023437500000000000"
     soln = -4.099527305346944316966073529329150915145874023437500000000000;
     ans = quantile_t(p=p, df=df);

     eps = 1e-50;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));

     // q <- qt(0.01, 3.45, lower.tail=F)
     // sprintf("%.60f", q)
     // "4.099527305346944316966073529329150915145874023437500000000000"
     soln = 4.099527305346944316966073529329150915145874023437500000000000;
     ans = quantile_t(p=p, df=df, false);

     eps = 1e-50;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));

}

//=========================================================================//


TEST_CASE("cdf chisq values"){

    // x <- pchisq(0.54, df=2.0, lower.tail=T, log.p=F)
    // sprintf("%.60f", x)
    // [1] "0.236620505663146823982501132377365138381719589233398437500000"

    // x <- pchisq(0.54, df=2.0, lower.tail=F, log.p=F)
    // sprintf("%.60f", x)
    // [1] "0.763379494336853148261923251993721351027488708496093750000000"
    
    // x <- pchisq(0.54, df=2.0, lower.tail=T, log.p=T)
    // sprintf("%.60f", x)
    // [1] "-1.441297663132672601804529222135897725820541381835937500000000"

    // x <- pchisq(0.54, df=2.0, lower.tail=F, log.p=T)
    // sprintf("%.60f", x)
    // [1] "-0.270000000000000017763568394002504646778106689453125000000000"


    double x = 0.54;
    double df = 2.0;
    double r_p1 = 0.236620505663146823982501132377365138381719589233398437500000;
    double r_p1_up = 0.763379494336853148261923251993721351027488708496093750000000;
    double r_p1_log = -1.441297663132672601804529222135897725820541381835937500000000;
    double r_p1_up_log= -0.270000000000000017763568394002504646778106689453125000000000;

    double e_cpp_p1 = cdf_chisq(x, df, true);
    double e_cpp_p1_log = cdf_chisq_log(x, df, true);
    double e_cpp_p1_up = cdf_chisq(x, df, false);
    double e_cpp_p1_up_log = cdf_chisq_log(x, df, false);


    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}




TEST_CASE("quantile chisq values"){
    double p = 0.95;
    double df = 12.3;
    // q <- qchisq(0.95, df=12.3)
    // sprintf("%.60f", q)
    // "21.428342879876034032804454909637570381164550781250000000000000"
    double soln = 21.428342879876034032804454909637570381164550781250000000000000;
    double ans = quantile_chisq(p=p, df=df);

    double eps = 1e-50;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));


     p = 0.01;
     df = 3.45;
     // q <- qchisq(0.01, df=3.45)
     // sprintf("%.60f", q)
     // "0.186609680758888452078991804228280670940876007080078125000000"
     soln = 0.186609680758888452078991804228280670940876007080078125000000;
     ans = quantile_chisq(p=p, df=df);

     eps = 1e-50;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));

     // q <- qchisq(0.01, 3.45, lower.tail=F)
     // sprintf("%.60f", q)
     // "12.232943156648634186467461404390633106231689453125000000000000"
     soln = 12.232943156648634186467461404390633106231689453125000000000000; 
     ans = quantile_chisq(p=p, df=df, false);

     eps = 1e-15;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));
}


//=========================================================================//

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

    double e_cpp_p1 = cdf_normal(x, mu, sigma);
    double e_cpp_p1_log = cdf_normal_log(x, mu, sigma,  true);
    double e_cpp_p1_up = cdf_normal(x, mu, sigma, false);
    double e_cpp_p1_up_log = cdf_normal_log(x, mu, sigma, false);

    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}

//=========================================================================//


TEST_CASE("quantile norm values"){
    double p = 0.95;
    double mu = 0;
    double sigma = 1;
    // R code
    // > q <- qnorm(0.95, 0, 1)
    // > sprintf("%.60f", q)
    // [1] "1.644853626951471525785564153920859098434448242187500000000000"
    double soln = 1.644853626951471525785564153920859098434448242187500000000000;
    double ans = quantile_normal(p=p, mu=mu, sigma=sigma);

    double eps = 1e-50;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));


     p = 0.01;
     mu = 0;
     sigma = 1;
     // > q <- qnorm(0.01, 0, 1)
     // > sprintf("%.60f", q)
     // [1] "-2.326347874040840757459136511897668242454528808593750000000000"
     soln = -2.326347874040840757459136511897668242454528808593750000000000;
     ans = quantile_normal(p, mu, sigma);

     eps = 1e-50;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));

}

//=========================================================================//


TEST_CASE("cdf gamma values"){

    // p <- pgamma(0.12, shape=1.2, scale=2.3)
    // sprintf("%.60f", p)
    // "0.025499214669880032602122810203582048416137695312500000000000"
    //
    // p <- pgamma(0.12, shape=1.2, scale=2.3, lower.tail=F, log.p=F)
    // sprintf("%.60f", p)
    // "0.974500785330120078420179652312071993947029113769531250000000"
    //
    // p <- pgamma(0.12, shape=1.2, scale=2.3, lower.tail=T, log.p=T)
    // sprintf("%.60f", p)
    // "-3.669107624551608548557624089880846440792083740234375000000000"
    //
    // p <- pgamma(0.12, shape=1.2, scale=2.3, lower.tail=F, log.p=T)
    // sprintf("%.60f", p)
    // "-0.025829954154784746239714721127711527515202760696411132812500"


    double p = 0.12;
    double shape = 1.2;
    double scale = 2.3;
    double r_p1 = 0.025499214669880032602122810203582048416137695312500000000000;
    double r_p1_up = 0.974500785330120078420179652312071993947029113769531250000000;
    double r_p1_log = -3.669107624551608548557624089880846440792083740234375000000000;
    double r_p1_up_log= -0.025829954154784746239714721127711527515202760696411132812500;

    double e_cpp_p1 = cdf_gamma(p, shape, scale, true);
    double e_cpp_p1_log = cdf_gamma_log(p, shape, scale, true);
    double e_cpp_p1_up = cdf_gamma(p, shape, scale, false);
    double e_cpp_p1_up_log = cdf_gamma_log(p, shape, scale, false);


    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}



TEST_CASE("quantile gamma values"){
    double p = 0.95;
    double shape = 3.4;
    double scale = 5.6;
    // q <- qgamma(0.95, shape=3.4, scale=5.6)
    // sprintf("%.60f", q)
    // "38.570354345912015503472503041848540306091308593750000000000000"
    double soln = 38.570354345912015503472503041848540306091308593750000000000000;
    double ans = quantile_gamma(p, shape, scale);

    double eps = 1e-50;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));


     p = 0.01;
     shape = 2.3;
     scale = 7.8;
     // q <- qgamma(0.01, shape=2.3, scale=7.8)
     // sprintf("%.60f", q)
     // "1.729272651898179802287813799921423196792602539062500000000000"
     soln = 1.729272651898179802287813799921423196792602539062500000000000;
     ans = quantile_gamma(p, shape, scale);

     eps = 1e-15;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));

     // q <- qgamma(0.01, shape=2.3, scale=7.8, lower.tail=F)
     // sprintf("%.60f", q)
     // "56.059175083606355372012330917641520500183105468750000000000000"
     soln = 56.059175083606355372012330917641520500183105468750000000000000; 
     ans = quantile_gamma(p, shape, scale, false);

     eps = 1e-15;
     CHECK(ans == doctest::Approx(soln).epsilon(eps));
}

//=========================================================================//


TEST_CASE("cdf wilcoxon values"){
    double q = 73;
    double m = 10;
    double n = 12;
    // set.seed(1)
    // x <- rnorm(10)
    // y <- rnorm(12, mean=1)
    // r <- rank(c(x, y))
    // r
    // [1]  4  6  2 16  7  3  9 11 10  5 22 14  8  1 21 12 13 20 18 15 19 17
    // sum(r[1:10])
    // [1] 73
    // p <- pwilcox(73, 10, 12)
    // sprintf("%.60f", p)
    // "0.809514015396368247223790604039095342159271240234375000000000"
    //
    // p <- pwilcox(73, 10, 12, F)
    // sprintf("%.60f", p)
    // "0.190485984603631697265058164703077636659145355224609375000000"
    //
    // p <- pwilcox(73, 10, 12, T, T)
    // sprintf("%.60f", p)
    // "-0.211321192368411542306105843636032659560441970825195312500000"
    //
    // p <- pwilcox(73, 10, 12, F, T)
    // sprintf("%.60f", p)
    // "-1.658176658756386201432064808614086359739303588867187500000000"

    double r_p1 = 0.809514015396368247223790604039095342159271240234375000000000;
    double r_p1_up = 0.190485984603631697265058164703077636659145355224609375000000;
    double r_p1_log = -0.211321192368411542306105843636032659560441970825195312500000;
    double r_p1_up_log= -1.658176658756386201432064808614086359739303588867187500000000;

    double e_cpp_p1 = cdf_wilcoxon(q, m, n);
    double e_cpp_p1_log = cdf_wilcoxon_log(q, m, n);
    double e_cpp_p1_up = cdf_wilcoxon(q, m, n, false);
    double e_cpp_p1_up_log = cdf_wilcoxon_log(q, m, n, false);

    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}


//=========================================================================//

TEST_CASE("cdf  noncentralchisq values"){

    // > x <- pchisq(0.71, df=2.48, ncp=1.23, lower.tail=T, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.120198917857169643164105821142584318295121192932128906250000"

    // > x <- pchisq(0.71, df=2.48, ncp=1.23, lower.tail=F, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.879801082142830370713681986671872437000274658203125000000000"
    
    // > x <- pchisq(0.71, df=2.48, ncp=1.23, lower.tail=T, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-2.118607259773721995799178330344147980213165283203125000000000"

    // > x <- pchisq(0.71, df=2.48, ncp=1.23, lower.tail=F, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-0.128059440080969755282325195366865955293178558349609375000000"


    double x = 0.71;
    double df = 2.48;
    double ncp = 1.23;
    double r_p1 = 0.120198917857169643164105821142584318295121192932128906250000;
    double r_p1_up = 0.879801082142830370713681986671872437000274658203125000000000;
    double r_p1_log = -2.118607259773721995799178330344147980213165283203125000000000;
    double r_p1_up_log= -0.128059440080969755282325195366865955293178558349609375000000;

    double e_cpp_p1 = cdf_noncentral_chisq(x, df, ncp, true);
    double e_cpp_p1_log = cdf_noncentral_chisq_log(x, df, ncp, true);
    double e_cpp_p1_up = cdf_noncentral_chisq(x, df, ncp, false);
    double e_cpp_p1_up_log = cdf_noncentral_chisq_log(x, df, ncp, false);


    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}




TEST_CASE("quantile chisq values"){
    double p = 0.95;
    double df = 12.3;
    double ncp = 3.69;
    // > q <- qchisq(p=0.95, df=12.3, ncp=3.69)
    // > sprintf("%.60f", q)
    // [1] "27.483655786335589255031663924455642700195312500000000000000000"
    double soln = 27.483655786335589255031663924455642700195312500000000000000000;
    double ans = quantile_noncentral_chisq(p=p, df=df, ncp=ncp);

    double eps = 1e-50;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));


    p = 0.01;
    // > q <- qchisq(p=0.01, df=12.3, ncp=3.69)
    // > sprintf("%.60f", q)
    // [1] "4.970213549946586262251457810634747147560119628906250000000000"
    soln = 4.970213549946586262251457810634747147560119628906250000000000;
    ans = quantile_noncentral_chisq(p=p, df=df, ncp=ncp);

    eps = 1e-50;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));

    // > q <- qchisq(p=0.01, df=12.3, ncp=3.69, lower.tail=F)
    // > sprintf("%.60f", q)
    // [1] "33.932873249693969341933552641421556472778320312500000000000000"
    soln = 33.932873249693969341933552641421556472778320312500000000000000; 
    ans = quantile_noncentral_chisq(p=p, df=df, ncp=ncp, false);

    eps = 1e-15;
    CHECK(ans == doctest::Approx(soln).epsilon(eps));
}


//=========================================================================//


TEST_CASE("cdf beta distribution"){

    // > x <- pbeta(0.7, 1.2, 3.4, lower.tail=T, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.977257254806371933320008338341722264885902404785156250000000"

    // > x <- pbeta(0.7, 1.2, 3.4, lower.tail=T, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-0.023005350641185102339436596707855642307549715042114257812500"

    // > x <- pbeta(0.7, 1.2, 3.4, lower.tail=F, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.022742745193628080557779469472734490409493446350097656250000"

    // > x <- pbeta(0.7, 1.2, 3.4, lower.tail=F, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-3.783509077030007716757609159685671329498291015625000000000000"


    double x = 0.7;
    double a = 1.2;
    double b = 3.4;
    double r_p1 = 0.977257254806371933320008338341722264885902404785156250000000;
    double r_p1_log = -0.023005350641185102339436596707855642307549715042114257812500;
    double r_p1_up = 0.022742745193628080557779469472734490409493446350097656250000;
    double r_p1_up_log= -3.783509077030007716757609159685671329498291015625000000000000;


    double e_cpp_p1 = cdf_beta(x, a, b, true);
    double e_cpp_p1_log = cdf_beta_log(x, a, b, true);
    double e_cpp_p1_up = cdf_beta(x, a, b, false);
    double e_cpp_p1_up_log = cdf_beta_log(x, a, b, false);
    
    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}//testing beta


//=========================================================================//

TEST_CASE("cdf Poisson distribution"){

    // > x <- ppois(q=8, lambda=7.89, lower.tail=T, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.607897963393453144576028535084333270788192749023437500000000"

    // > x <- ppois(q=8, lambda=7.89, lower.tail=T, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-0.497748234465917704927306886020232923328876495361328125000000"

    // > x <- ppois(q=8, lambda=7.89, lower.tail=F, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.392102036606546799912820233657839708030223846435546875000000"

    // > x <- ppois(q=8, lambda=7.89, lower.tail=F, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-0.936233175597501743325778988946694880723953247070312500000000"


    double x = 8;
    double lambda = 7.89;
    double r_p1 = 0.607897963393453144576028535084333270788192749023437500000000;
    double r_p1_log = -0.497748234465917704927306886020232923328876495361328125000000;
    double r_p1_up = 0.392102036606546799912820233657839708030223846435546875000000;
    double r_p1_up_log= -0.936233175597501743325778988946694880723953247070312500000000;


    double e_cpp_p1 = cdf_poisson(x, lambda, true);
    double e_cpp_p1_log = cdf_poisson_log(x, lambda, true);
    double e_cpp_p1_up = cdf_poisson(x, lambda, false);
    double e_cpp_p1_up_log = cdf_poisson_log(x, lambda, false);
    
    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}//testing Poisson 

//=========================================================================//

TEST_CASE("cdf Binomial distribution"){

    // > x <- pbinom(q=8, size=20, prob=0.7, lower.tail=T, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.005138161535121414658089378235672484152019023895263671875000"

    // > x <- pbinom(q=8, size=20, prob=0.7, lower.tail=T, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-5.271059941489036226869302481645718216896057128906250000000000"


    // > x <- pbinom(q=8, size=20, prob=0.7, lower.tail=F, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.994861838464878633914167949114926159381866455078125000000000"

    // > x <- pbinom(q=8, size=20, prob=0.7, lower.tail=F, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-0.005151407279097742521190017583876397111453115940093994140625"


    int x = 8;
    int n = 20;
    double p = 0.7;
    double r_p1 = 0.005138161535121414658089378235672484152019023895263671875000;
    double r_p1_log = -5.271059941489036226869302481645718216896057128906250000000000;
    double r_p1_up = 0.994861838464878633914167949114926159381866455078125000000000;
    double r_p1_up_log= -0.005151407279097742521190017583876397111453115940093994140625;


    double e_cpp_p1 = cdf_binomial(x, n, p, true);
    double e_cpp_p1_log = cdf_binomial_log(x, n, p, true);
    double e_cpp_p1_up = cdf_binomial(x, n, p, false);
    double e_cpp_p1_up_log = cdf_binomial_log(x, n, p, false);
    
    double eps = 1e-30;
    double eps_log = 1e-18;

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}//testing binomial 


//=========================================================================//


TEST_CASE("cdf exponential distribution"){

    // > x <- pexp(q=2.34, rate=1.23, lower.tail=T, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.943764103599861514659608019428560510277748107910156250000000"

    // > x <- pexp(q=2.34, rate=1.23, lower.tail=T, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-0.057879034318792439706147234801392187364399433135986328125000"

    // > x <- pexp(q=2.34, rate=1.23, lower.tail=F, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.056235896400138499218179788385896245017647743225097656250000"


    // > x <- pexp(q=2.34, rate=1.23, lower.tail=F, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-2.878200000000000091660012913052923977375030517578125000000000"

    double x = 2.34;
    double rate = 1.23;
    double r_p1 = 0.943764103599861514659608019428560510277748107910156250000000;
    double r_p1_log = -0.057879034318792439706147234801392187364399433135986328125000;
    double r_p1_up = 0.056235896400138499218179788385896245017647743225097656250000;
    double r_p1_up_log= -2.878200000000000091660012913052923977375030517578125000000000;


    double e_cpp_p1 = cdf_exp(x, rate, true);
    double e_cpp_p1_log = cdf_exp_log(x, rate, true);
    double e_cpp_p1_up = cdf_exp(x, rate, false);
    double e_cpp_p1_up_log = cdf_exp_log(x, rate, false);
    
    double eps = 1e-30;
    double eps_log = 1e-18;

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
} // exponential distribution

//=========================================================================//


TEST_CASE("cdf geometric distribution"){

    // > x <- pgeom(q=3, p=0.6, lower.tail=T, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.974400000000000043876013933186186477541923522949218750000000"

    // > x <- pgeom(q=3, p=0.6, lower.tail=T, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-0.025933382026504483985895888054074021056294441223144531250000"

    // > x <- pgeom(q=3, p=0.6, lower.tail=F, log.p=F)
    // > sprintf("%.60f", x)
    // [1] "0.025600000000000008165690346118026354815810918807983398437500"

    // > x <- pgeom(q=3.1, p=0.6, lower.tail=F, log.p=T)
    // > sprintf("%.60f", x)
    // [1] "-3.665162927496619982292713757487945258617401123046875000000000"

    double x = 3;
    double p = 0.6;
    double r_p1 = 0.974400000000000043876013933186186477541923522949218750000000;
    double r_p1_log = -0.025933382026504483985895888054074021056294441223144531250000;
    double r_p1_up = 0.025600000000000008165690346118026354815810918807983398437500;
    double r_p1_up_log= -3.665162927496619982292713757487945258617401123046875000000000;


    double e_cpp_p1 = cdf_geom(x, p, true);
    double e_cpp_p1_log = cdf_geom_log(x, p, true);
    double e_cpp_p1_up = cdf_geom(x, p, false);
    double e_cpp_p1_up_log = cdf_geom_log(x, p, false);
    
    double eps = 1e-30;
    double eps_log = 1e-17;

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
} // exponential distribution

//=========================================================================//


TEST_CASE("cdf Cauchy values"){
    double x = 3.45;
    double loc = 1.2;
    double scale = 1.5;

    // > p <- pcauchy(3.45, loc=1.2, scale=1.5, lower.tail=T, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.812832958189001253401784197194501757621765136718750000000000"
    
    // > p <- pcauchy(3.45, loc=1.2, scale=1.5, lower.tail=F, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.187167041810998802109367034063325263559818267822265625000000"
    
    // > p <- pcauchy(3.45, loc=1.2, scale=1.5, lower.tail=T, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-0.207229654027002746508046016060688998550176620483398437500000"
    
    // > p <- pcauchy(3.45, loc=1.2, scale=1.5, lower.tail=F, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-1.675753789140727478823578167066443711519241333007812500000000"


    double r_p1 = 0.812832958189001253401784197194501757621765136718750000000000;
    double r_p1_up = 0.187167041810998802109367034063325263559818267822265625000000;
    double r_p1_log = -0.207229654027002746508046016060688998550176620483398437500000;
    double r_p1_up_log= -1.675753789140727478823578167066443711519241333007812500000000;

    double e_cpp_p1 = cdf_cauchy(x, loc, scale);
    double e_cpp_p1_log = cdf_cauchy_log(x, loc, scale,  true);
    double e_cpp_p1_up = cdf_cauchy(x, loc, scale, false);
    double e_cpp_p1_up_log = cdf_cauchy_log(x, loc, scale, false);

    double eps = 1e-50;
    double eps_log = 1e-16;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
} // Cauchy distribution


//=========================================================================//
TEST_CASE("cdf weibull values"){

    // > p <- pweibull(0.69, shape=1.2, scale=2.3, lower.tail=T, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.210062085365011030901882804755587130784988403320312500000000"
    //
    // > p <- pweibull(0.69, shape=1.2, scale=2.3, lower.tail=F, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.789937914634988969098117195244412869215011596679687500000000"
    //
    // > p <- pweibull(0.69, shape=1.2, scale=2.3, lower.tail=T, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-1.560352147363629260468087522895075380802154541015625000000000"
    //
    // > p <- pweibull(0.69, shape=1.2, scale=2.3, lower.tail=F, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-0.235800925678986833533556932707142550498247146606445312500000"


    double p = 0.69;
    double shape = 1.2;
    double scale = 2.3;
    double r_p1 = 0.210062085365011030901882804755587130784988403320312500000000;
    double r_p1_up = 0.789937914634988969098117195244412869215011596679687500000000;
    double r_p1_log = -1.560352147363629260468087522895075380802154541015625000000000;
    double r_p1_up_log= -0.235800925678986833533556932707142550498247146606445312500000;

    double e_cpp_p1 = cdf_weibull(p, shape, scale, true);
    double e_cpp_p1_log = cdf_weibull_log(p, shape, scale, true);
    double e_cpp_p1_up = cdf_weibull(p, shape, scale, false);
    double e_cpp_p1_up_log = cdf_weibull_log(p, shape, scale, false);


    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}
//=========================================================================//


TEST_CASE("cdf hypergeometric"){

    // > p <- phyper(4, m=5, n=6, k=7, lower.tail=T, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.954545454545454585826291804551146924495697021484375000000000"
    //
    // > p <- phyper(4, m=5, n=6, k=7, lower.tail=F, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.045454545454545393357026483727167942561209201812744140625000"
    //
    // > p <- phyper(4, m=5, n=6, k=7, lower.tail=T, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-0.046520015634892809830436277707121917046606540679931640625000"
    //
    // > p <- phyper(4, m=5, n=6, k=7, lower.tail=F, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-3.091042453358316954847850865917280316352844238281250000000000"


    double x = 4;
    double m = 5;
    double n = 6; 
    double k = 7;
    double r_p1 = 0.954545454545454585826291804551146924495697021484375000000000;
    double r_p1_up = 0.045454545454545393357026483727167942561209201812744140625000;
    double r_p1_log = -0.046520015634892809830436277707121917046606540679931640625000;
    double r_p1_up_log= -3.091042453358316954847850865917280316352844238281250000000000;

    double e_cpp_p1 = cdf_hypergeometric(x, m, n, k, true);
    double e_cpp_p1_log = cdf_hypergeometric_log(x, m, n, k, true);
    double e_cpp_p1_up = cdf_hypergeometric(x, m, n, k, false);
    double e_cpp_p1_up_log = cdf_hypergeometric_log(x, m, n, k, false);


    double eps = 1e-16;
    double eps_log = 1e-17;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}
//=========================================================================//

TEST_CASE("cdf lognormal"){

    // > p <- plnorm(1.23, meanlog=0, sdlog=1,lower.tail=T, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.582000603689423301467797955410787835717201232910156250000000"
    //
    // > p <- plnorm(1.23, meanlog=0, sdlog=1,lower.tail=F, log.p=F)
    // > sprintf("%.60f", p)
    // [1] "0.417999396310576698532202044589212164282798767089843750000000"
    //
    // > p <- plnorm(1.23, meanlog=0, sdlog=1,lower.tail=T, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-0.541283793984186800685165508184581995010375976562500000000000"
    //
    // > p <- plnorm(1.23, meanlog=0, sdlog=1,lower.tail=F, log.p=T)
    // > sprintf("%.60f", p)
    // [1] "-0.872275290691493787598176368192071095108985900878906250000000"

    double x = 1.23;
    double meanlog = 0;
    double sdlog = 1; 
    double r_p1 = 0.582000603689423301467797955410787835717201232910156250000000;
    double r_p1_up = 0.417999396310576698532202044589212164282798767089843750000000;
    double r_p1_log = -0.541283793984186800685165508184581995010375976562500000000000;
    double r_p1_up_log= -0.872275290691493787598176368192071095108985900878906250000000;

    double e_cpp_p1 = cdf_lognormal(x, meanlog, sdlog, true);
    double e_cpp_p1_log = cdf_lognormal_log(x, meanlog, sdlog, true);
    double e_cpp_p1_up = cdf_lognormal(x, meanlog, sdlog, false);
    double e_cpp_p1_up_log = cdf_lognormal_log(x, meanlog, sdlog, false);


    double eps = 1e-16;
    double eps_log = 1e-17;

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}
//=========================================================================//
