#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

//include others
#include <iostream>
#include "../cdf_base.cpp"

//global boolean flag
bool SHOW_STRING = false;
// bool SHOW_STRING = true;
//

TEST_CASE("Example showing tolerance"){
    REQUIRE( (2.1 - 1e-10) == doctest::Approx(2.1));

    // allow for a 1% error
    REQUIRE(22.0/7 == doctest::Approx(3.141).epsilon(0.01)); 
}

TEST_CASE("First constants"){
    //CHECK(2 ==2*1);
    CHECK(ML_POSINF == ML_POSINF);
    if (SHOW_STRING)
        std::cout << " 1. ML_POSINF: " << ML_POSINF << std::endl;
    
    CHECK(ML_NEGINF == ML_NEGINF);
    if (SHOW_STRING)
        std::cout << " 2. ML_NEGINF: " << ML_NEGINF << std::endl;

    CHECK(ISNAN(ML_NAN)== true);
    if (SHOW_STRING)
        std::cout << " 3. NAN: " << ML_NAN << std::endl;

}



TEST_CASE("R_D_ functions"){
    //need these two defined for the macro
    bool log_p = true;
    bool lower_tail = false;

    CHECK(R_D__0 == R_D__0);
    if (SHOW_STRING)
        std::cout << " 4. R_D__0: " << R_D__0 << std::endl;

    CHECK(R_D__1 == R_D__1);
    if (SHOW_STRING)
        std::cout << " 5. R_D__1: " << R_D__1 << std::endl;

    CHECK(R_DT_0 == R_DT_0);
    if (SHOW_STRING)
        std::cout << " 6. R_DT_0: " << R_DT_0 << std::endl;

    CHECK(R_DT_1 == R_DT_1);
    if (SHOW_STRING)
        std::cout << " 7. R_DT_1: " << R_DT_1 << std::endl;
}


TEST_CASE("Second contants"){
    CHECK(ME_PRECISION == ME_PRECISION);
    if (SHOW_STRING)
        std::cout << " 8. ME_PRECISION: " << ME_PRECISION << std::endl;

    CHECK(ME_UNDERFLOW == ME_UNDERFLOW);
    if (SHOW_STRING)
        std::cout << " 9. ME_UNDERFLOW: " << ME_UNDERFLOW << std::endl;

    CHECK(ME_RANGE == ME_RANGE);
    if (SHOW_STRING)
        std::cout << "10. ME_RANGE: " << ME_RANGE << std::endl;

    CHECK(ME_DOMAIN == ME_DOMAIN);
    if (SHOW_STRING)
        std::cout << "11. ME_DOMAIN: " << ME_DOMAIN << std::endl;

    CHECK(ME_NOCONV == ME_NOCONV);
    if (SHOW_STRING)
        std::cout << "12. ME_NOCONV: " << ME_NOCONV << std::endl;

    CHECK(scalefactor == scalefactor);
    if (SHOW_STRING)
        std::cout << "13. scalefactor: " << scalefactor << std::endl;

}


TEST_CASE("Checking DBL_ constants"){
    //already in float.h
    CHECK(DBL_EPSILON == DBL_EPSILON);
    if (SHOW_STRING)
        std::cout << "14. DBL_EPSILON: " << DBL_EPSILON << std::endl;

    CHECK(DBL_MAX == DBL_MAX);
    if (SHOW_STRING)
        std::cout << "15. DBL_MAX: " << DBL_MAX << std::endl;

    CHECK(DBL_MIN == DBL_MIN);
    if (SHOW_STRING)
        std::cout << "16. DBL_MIN: " << DBL_MIN << std::endl;

}


TEST_CASE("Checking M_ constantts"){
    //already in math.h
    //but have #indef in cdf_base.cpp

    CHECK(M_PI == M_PI);
    if (SHOW_STRING)
        std::cout << "17. M_PI: " << M_PI << std::endl;

    CHECK(M_1_SQRT_2PI == M_1_SQRT_2PI);
    if (SHOW_STRING)
        std::cout << "18. M_1_SQRT_2PI: " << M_1_SQRT_2PI << std::endl;

//     CHECK(M_LN_SQRT_PI == M_LN_SQRT_PI);
//     if (SHOW_STRING)
//         std::cout << "18. M_LN_SQRT_PI: " << M_LN_SQRT_PI << std::endl;

    CHECK(M_LN_SQRT_2PI == M_LN_SQRT_2PI);
    if (SHOW_STRING)
        std::cout << "19. M_LN_SQRT_2PI: " << M_LN_SQRT_2PI << std::endl;

    CHECK(M_LN_SQRT_PId2 == M_LN_SQRT_PId2);
    if (SHOW_STRING)
        std::cout << "20. M_LN_SQRT_PId2: " << M_LN_SQRT_PId2 << std::endl;

    CHECK(M_LN2 == M_LN2);
    if (SHOW_STRING)
        std::cout << "21. M_LN2: " << M_LN2 << std::endl;

    CHECK(M_LOG10_2 == M_LOG10_2);
    if (SHOW_STRING)
        std::cout << "22. M_LOG10_2: " << M_LOG10_2 << std::endl;

    CHECK(M_SQRT_32 == M_SQRT_32);
    if (SHOW_STRING)
        std::cout << "23. M_SQRT_32: " << M_SQRT_32 << std::endl;
    
    CHECK(M_2PI == M_2PI);
    if (SHOW_STRING)
        std::cout << "24. M_2PI: " << M_2PI << std::endl;
}

TEST_CASE("Native to CPP"){
    CHECK(CHAR_BIT == CHAR_BIT);
    if (SHOW_STRING)
        std::cout << "25. CHAR_BIT: " << CHAR_BIT << std::endl;

    CHECK(INT_MAX == INT_MAX);
    if (SHOW_STRING)
        std::cout << "26. INT_MAX: " << INT_MAX << std::endl;
}



TEST_CASE("FLT_ constants"){
    //Included in math.h
    //actually included in float.h?
    CHECK(FLT_RADIX == FLT_RADIX);
    if (SHOW_STRING)
        std::cout << "27. FLT_RADIX: " << FLT_RADIX << std::endl;

    CHECK(FLT_MANT_DIG == FLT_MANT_DIG);
    if (SHOW_STRING)
        std::cout << "28. FLT_MANT_DIG: " << FLT_MANT_DIG << std::endl;

    CHECK(FLT_MIN_EXP == FLT_MIN_EXP);
    if (SHOW_STRING)
        std::cout << "29. FLT_MIN_EXP: " << FLT_MIN_EXP << std::endl;

    CHECK(FLT_MAX_EXP == FLT_MAX_EXP);
    if (SHOW_STRING)
        std::cout << "30. FLT_MAX_EXP: " << FLT_MAX_EXP << std::endl;
}



TEST_CASE("More DBL_ constants"){
    CHECK(DBL_MANT_DIG == DBL_MANT_DIG);
    if (SHOW_STRING)
        std::cout << "31. DBL_MANT_DIG: " << DBL_MANT_DIG << std::endl;

    CHECK(DBL_MIN_EXP == DBL_MIN_EXP);
    if (SHOW_STRING)
        std::cout << "32. DBL_MIN_EXP: " << DBL_MIN_EXP << std::endl;

    CHECK(DBL_MAX_EXP == DBL_MAX_EXP);
    if (SHOW_STRING)
        std::cout << "33. DBL_MAX_EXP: " << DBL_MAX_EXP << std::endl;
}


TEST_CASE("Another constant"){
    CHECK(M_cutoff == M_cutoff);
    if (SHOW_STRING)
        std::cout << "34. M_cutoff: " << M_cutoff << std::endl;
}

//From now on, no more printing to screen, just checking functions


//the problem is that R_P_bounds_01 is a macro with a return,
//so this function forces the return
double check_R_P_bounds01(double x){
    //need the lower_tail and log_p booleans defined for the macro
    bool lower_tail = true; 
    bool log_p = true;
    R_P_bounds_01(x, 0., ML_POSINF);
    return 0;
}


TEST_CASE("checking UTILS function R_P_bounds_01"){
    //WEIRDLY need the lower_tail and log_p booleans here too
    //to define R_DT_0, etc
    bool lower_tail = true; 
    bool log_p = true;
    CHECK(check_R_P_bounds01(-1) == R_DT_0);
    CHECK(check_R_P_bounds01(3) == 0);
    CHECK(check_R_P_bounds01(ML_POSINF) == R_DT_1);
}


TEST_CASE("checking more UTILS functions"){
    CHECK(R_FINITE(1)==true);
    CHECK(R_FINITE(ML_POSINF)==false);
    CHECK(R_FINITE(ML_NAN)==false);
    CHECK(R_FINITE(ML_NEGINF)==false);

    CHECK(R_forceint(1.2)==1);

    bool log_p = true;
    bool lower_tail = true;
    //give_log is an alias for log_p
    //bool give_log = true;
    double p = 0.01;
    CHECK(R_D_Cval(p) == (1-p));
    CHECK(R_D_exp(p)==p);
    double f = 0.05;
    double x = 0.2;
    double ans = -0.5*log(f) + x;
    CHECK(R_D_fexp(f, x) == ans);
}


TEST_CASE("checking Rboolean"){
    Rboolean doswap=true;
    Rboolean donotswap=false;
//     Rboolean doswap=TRUE;
//     Rboolean donotswap=FALSE;

    //the problem was that this would not work with
    //the original enum definition of Rboolean
    //it would not cast from bool to Rboolean
    double x = 1.1;
    doswap = (x > 0.5);

    CHECK(doswap==true);
    CHECK(donotswap==false);
}


TEST_CASE("Native C++ in math.h"){
    CHECK(fmod(-2.2, 2.0)==doctest::Approx(-0.2) );
    CHECK(fmod(4.3, 2.0)==doctest::Approx(0.3) );
    //CHECK(fmod(4.3, 2.0)==0.3);

    CHECK(sin(0.5)==sin(0.5));
    CHECK(log(0.5)==log(0.5));
    CHECK(exp(0.5)==exp(0.5));
    CHECK(pow(0.5, 0.5)==pow(0.5, 0.5));
    CHECK(sqrt(4.0)==2.0);

    CHECK(fabs(-0.5)==0.5);
    CHECK(ldexp(4.0, -2) == 1.0);

    CHECK(lgamma(3)==lgamma(3));

    //checking <<
    int j = 1;
    CHECK((j<<1)==2);
    CHECK(j==1);
//     int y = (j<<1);
//     std::cout << j << " " << y << std::endl;

    CHECK(round(0.8)==1);
    CHECK(nearbyint(1.4)==1);

    CHECK(isfinite(1)==true);
    CHECK(isfinite(ML_NAN)==false);
}

TEST_CASE("Warnings"){
    double y = 0.0;
    MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
}

TEST_CASE("More elementary constants"){
    CHECK(Rf_d1mach(1)==DBL_MIN);

    //d1mach uses a pointer
    int j = 2;
    CHECK(d1mach(&j)==DBL_MAX);

    double x = 0.5;
    const double a = 0.1;
    const int n = 2;

    CHECK(chebyshev_eval(x, &a, n)==chebyshev_eval(x, &a, n));

    //now will be NAN
    x = 1.5;
    CHECK( ISNAN(chebyshev_eval(x, &a, n) ) );

// int chebyshev_init(double *dos, int nos, double eta){
    
    double dos = 0.5;
    int nos = 1;
    double eta = 0.01;
    CHECK( chebyshev_init(&dos, nos, eta) == chebyshev_init(&dos, nos, eta));

    //log1p
    CHECK(log1p(x) == log1p(x));

    //lgammacor - works with big numbers
    double y = 20;
    CHECK(lgammacor(y) == lgammacor(y));


    //sinpi
    CHECK(sinpi(y) == sinpi(y));

    //fmax2
    CHECK(fmax2(1.0, 2.0) == 2.0);

    //gammalims
    double xmax4 = 1.0;
    double xmin1 = 0.1;
    gammalims(&xmin1, &xmax4);
    CHECK(xmin1==doctest::Approx(-170.567).epsilon(0.01) );
    CHECK(xmax4==doctest::Approx(171.614).epsilon(0.01) );

    //lgammaf used in math.h, and used in strilerr
    CHECK(lgammaf(0.5) == lgammaf(0.5));

    //stirlerr
    CHECK(stirlerr(0.5)==stirlerr(0.5));

    //gammafn
    CHECK(gammafn(0.7) == gammafn(0.7));

    //lgammafn
    CHECK(lgammafn(0.5) == lgammafn(0.5));

    //comparing lgammaf to lgammafn
    CHECK(lgammafn(0.5) == doctest::Approx(lgammaf(0.5)) );

    //lbeta
    CHECK(lbeta(0.5, 0.6) == lbeta(0.5, 0.6));

    //bd0
    CHECK( bd0(0.3, 0.5) == bd0(0.3, 0.5) );

    //dpois_raw
    CHECK( dpois_raw(0.5, 0.6, true) ==  dpois_raw(0.5, 0.6, true) );

// dpois_wrap (double x_plus_1, double lambda, int give_log)

    //dpois_wrap
    CHECK(dpois_wrap(0.7, 0.8, false) == dpois_wrap(0.7, 0.8, false));

    //dnorm4 and dnorm
    x = 0.7;
    double mu = 0.1;
    double sigma = 1.1;
    CHECK(dnorm4(x, mu, sigma, true) == dnorm4(x, mu, sigma, true));
    CHECK(dnorm(x, mu, sigma, true) == dnorm(x, mu, sigma, true));

    //dpnorm
    double lp = 0.8;
    CHECK( dpnorm(x, false, lp) == dpnorm(x, false, lp));
    //std::cout <<  dpnorm(x, false, lp) << std::endl;

    //logcf
    x = 0.02;
    double i = 3;
    double d = 2;
    double eps = 1e-14;
    CHECK(logcf(x, i , d, eps) == logcf(x, i , d, eps));

    //log1pmx
    CHECK(log1pmx(x) == log1pmx(x));

    //expm1
    double z = 0.7;
    CHECK(expm1(z) == expm1(z));
    z = -0.6;
    CHECK(expm1(z) == expm1(z));
    //std::cout << expm1(z) << std::endl;
    

    //lgamma1p
    CHECK(lgamma1p(x) == lgamma1p(x));
}

TEST_CASE("Now another round of functions, building to ppois_asymp"){
    double x = -0.6;
    CHECK(R_Log1_Exp(x) == R_Log1_Exp(x));
    x = 0.7;
    CHECK(ISNAN(R_Log1_Exp(x)));

    //pgamma_smallx
    double alph = 0.7;
    CHECK(pgamma_smallx(x, alph, true, true));

    //pd_upper_series
    double y = 0.4;
    CHECK(pd_upper_series(x, y, true) == pd_upper_series(x, y, true));

    //pd_lower_cf
    double d = 0.9;
    CHECK(pd_lower_cf(y, d) == pd_lower_cf(y, d));

    //pd_lower_series
    double lambda = 0.95;
    CHECK(pd_lower_series(lambda, y) == pd_lower_series(lambda, y));


    //pnorm_both
    x = 0.5;
    double cum = 2.0;
    double ccum = 3.2;
    int i_tail = 2;
    //changing values in place
    pnorm_both(x, &cum, &ccum, i_tail, true);

    //checking pnorm and pnorm5 at the same time
    x = 0.5;
    double mu = 1.2;
    double sigma = 2.3;
    CHECK(pnorm(x, mu, sigma, true, true)==pnorm5(x, mu, sigma, true, true));

    //ppois_asymp
    CHECK(ppois_asymp(x, 0.7, true, true) == ppois_asymp(x, 0.7, true, true));
}


TEST_CASE("Now for pgamma"){

    //pgamma_raw
    double x = 0.7;
    double a = 0.5;
    int lt = 0;
    int log_p = 0;
    CHECK(pgamma_raw(x, a, lt, log_p)==pgamma_raw(x, a, lt, log_p));

    //pgamma
    double s = 1.5;
    CHECK(pgamma(x, a, s, lt, log_p)==pgamma(x, a, s, lt, log_p));

    //pchisq
    double df = 3.2;
    CHECK(pchisq(x, df, lt, log_p)==pchisq(x, df, lt, log_p));
}



TEST_CASE("Elementary functions for pt"){

    //Rf_i1mach
    int i = 2;
    CHECK(Rf_i1mach(i) == Rf_i1mach(i));

    //exparg
    int j = 0;
    CHECK(exparg(j)==exparg(j));

    //fpser
    double a = 0.5;
    double b = 0.6;
    double x = 0.7;
    double eps = 0.001;
    CHECK(fpser(a, b, x, eps, true) == fpser(a, b, x, eps, true));

    //psi
    CHECK(psi(x)==psi(x));
    
    //apser
    CHECK(apser(a, b, x, eps) == apser(a, b, x, eps));
}



TEST_CASE("gam elementary functions"){
    //gam1
    double a = 0.8;
    CHECK(gam1(a) == gam1(a));

    //gamln1
    CHECK(gamln1(a) == gamln1(a));

    //gamln
    CHECK(gamln(a) == gamln(a));

    //alnrel
    CHECK(alnrel(a) == alnrel(a));

    //algdiv
    double b = 0.3;
    CHECK(algdiv(a, b) == algdiv(a, b));

    //bcorr
    CHECK(bcorr(a, b)==bcorr(a, b));

    //gsumln
    CHECK(gsumln(a, b)==gsumln(a, b));

    //betaln
    CHECK(betaln(a, b)==betaln(a, b));

    //bpser
    double x = 0.5;
    double eps = 0.001;
    CHECK(bpser(a, b, x, eps, true) == bpser(a, b, x, eps, true));
}


TEST_CASE("More elem functions"){

    //esum
    int mu = 3;
    double x = 0.5;
    CHECK(esum(mu, x, true) == esum(mu, x, true));

    //rlog1
    CHECK(rlog1(x)==rlog1(x));

    //brcmp1
    double a = 0.9;
    double b = 0.2;
    double y = 0.2;
    CHECK(brcmp1(mu, a, b, x, y, true)==brcmp1(mu, a, b, x, y, true));

    //rexpm1
    CHECK(rexpm1(x)==rexpm1(x));

    //erf__
    CHECK(erf__(x)==erf__(x));

    //erfc1
    int i = 1;
    CHECK(erfc1(i, x)==erfc1(i, x));

    //grat_r
    double lr = 0.1;
    double eps = 0.001;
    CHECK(grat_r(a, x, lr, eps)==grat_r(a, x, lr, eps));

    //bup
    int n = 2;
    CHECK(bup(a, b, x, y, n, eps, true)==bup(a, b, x, y, n, eps, true));

    //logspace_add
    CHECK(logspace_add(x, y)==logspace_add(x, y));


    //bgrat, no check, it is void
    double w = 0.3;
    int ierr = 2;
    bgrat(a, b, x, y, &w, eps, &ierr, true);

    //basym
    double la = 0.4;
    CHECK(basym(a, b, la, eps, true)==basym(a, b, la, eps, true));

    //brcomp
    CHECK(brcomp(a, b, x, y, true)==brcomp(a, b, x, y, true));

    //bfrac
    CHECK(bfrac(a, b, x, y, la, eps, true)==bfrac(a, b, x, y, la, eps, true));

    //bratio
    double w1 = 0.35;
    bratio(a, b, x, y, &w, &w1, &ierr, true);

}


TEST_CASE("pbeta functions"){

    //pbeta_raw
    double x = 0.5;
    double a = 0.3;
    double b = 0.7;
    CHECK(pbeta_raw(x, a, b, true, true)==pbeta_raw(x, a, b, true, true));

    //pbeta
    CHECK(pbeta(x, a, b, true, true)==pbeta(x, a, b, true, true));
}


TEST_CASE("pt function"){

    double x = 0.5;
    double n = 2.2;

    //pt
    CHECK(pt(x, n, true, true)==pt(x, n, true, true));
}



TEST_CASE("pnorm values"){
    //from R
    //pnorm(0.5, lower_tail=true, log.p=false)
    //pnorm(0.5, lower_tail=true, log.p=true)
    double r_p1 = 0.69146246127401300718418042379198595881462097167969;
    double r_p2 = 0.691462461;
    double r_p1_log = -0.36894641528865651514124124332738574594259262084961;
    // x <- pnorm(0.5, 0.0, 1.0, F, F)
    // sprintf("%.60f", x)
    //  [1] "0.308537538725986937304668344950187020003795623779296875000000"
    //
    //x <- pnorm(0.5, 0.0, 1.0, F, T)
    // sprintf("%.60f", x)
    // [1] "-1.175911761593618543031425360823050141334533691406250000000000"

    double r_p1_up = 0.308537538725986937304668344950187020003795623779296875000000;
    double r_p1_up_log = -1.175911761593618543031425360823050141334533691406250000000000;


    double q1 = 0.5;
    double mu = 0.0;
    double sigma = 1.0;
    double cpp_p1 = pnorm(q1, mu, sigma, true, false);
    double cpp_p1_log = pnorm(q1, mu, sigma, true, true);
    double cpp_p1_up = pnorm(q1, mu, sigma, false, false);
    double cpp_p1_up_log = pnorm(q1, mu, sigma, false, true);

    //want to check to 30 places
    double eps = 1e-50;
    double eps_log = 1e-50;
    CHECK(cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));

    //should fail
    CHECK(r_p2 != doctest::Approx(r_p1).epsilon(eps));

//     int prec = 100;
//     std::cout << std::setprecision(prec) << r_p1 << std::endl;
//     std::cout << std::setprecision(prec) << r_p2 << std::endl;
//     std::cout << std::setprecision(prec) << cpp_p1 << std::endl;
//     std::cout << std::setprecision(prec) << r_p1_log << std::endl;
//     std::cout << std::setprecision(prec) << cpp_p1_log << std::endl;

}


TEST_CASE("pgamma values"){
    // x <- pgamma(0.7, shape=1.2, scale=2.3, lower.tail=T, log.p=F)
    // sprintf("%.60f", x)
    // [1] "0.185100830776064373406342156158643774688243865966796875000000"
    //
    // x <- pgamma(0.7, shape=1.2, scale=2.3, lower.tail=T, log.p=T)
    // sprintf("%.60f", x)
    // [1] "-1.686854571157412108206585799052845686674118041992187500000000"
    //
    // x <- pgamma(0.7, shape=1.2, scale=2.3, lower.tail=F, log.p=F)
    // sprintf("%.60f", x)
    // [1] "0.814899169223935571082506612583529204130172729492187500000000"
    //
    // x <- pgamma(0.7, shape=1.2, scale=2.3, lower.tail=F, log.p=T)
    // sprintf("%.60f", x)
    // [1] "-0.204690892138706725944530262495391070842742919921875000000000"

    double x = 0.7;
    double shape = 1.2;
    double scale = 2.3;
    double r_p1 = 0.185100830776064373406342156158643774688243865966796875000000;
    double r_p1_log = -1.686854571157412108206585799052845686674118041992187500000000;
    double r_p1_up = 0.814899169223935571082506612583529204130172729492187500000000;
    double r_p1_up_log = -0.204690892138706725944530262495391070842742919921875000000000;

    double cpp_p1 = pgamma(x, shape, scale, true, false);
    double cpp_p1_log = pgamma(x, shape, scale, true, true);
    double cpp_p1_up = pgamma(x, shape, scale, false, false);
    double cpp_p1_up_log = pgamma(x, shape, scale, false, true);

    //will fail now
    double eps = 1e-50;
    double eps_log = 1e-50;
    CHECK(cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));



//     int prec = 100;
//     std::cout << std::setprecision(prec) << r_p1 << std::endl;
//     std::cout << std::setprecision(prec) << cpp_p1 << std::endl;
//     std::cout << std::setprecision(prec) << r_p1_log << std::endl;
//     std::cout << std::setprecision(prec) << cpp_p1_log << std::endl;
} //end of testing pgamma values


TEST_CASE("pt values"){

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
    double r_p1_up_log = -1.279427536176580559867943520657718181610107421875000000000000;



    double cpp_p1 = pt(x, df, true, false);
    double cpp_p1_log = pt(x, df, true, true);
    double cpp_p1_up = pt(x, df, false, false);
    double cpp_p1_up_log = pt(x, df, false, true);

    double e_cpp_p1 = cdf_t(x, df, true);
    double e_cpp_p1_log = cdf_t_log(x, df, true);
    double e_cpp_p1_up = cdf_t(x, df, false);
    double e_cpp_p1_up_log = cdf_t_log(x, df, false);
//     int prec = 100;
//     std::cout << std::setprecision(prec) << r_p1 << std::endl;
//     std::cout << std::setprecision(prec) << cpp_p1 << std::endl;
    
    double eps = 1e-50;
    double eps_log = 1e-50;
    CHECK(cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}//testing pt



TEST_CASE("pchisq values"){

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
    double r_p1_up_log = -0.270000000000000017763568394002504646778106689453125000000000;

    double cpp_p1 = pchisq(x, df, true, false);
    double cpp_p1_log = pchisq(x, df, true, true);
    double cpp_p1_up = pchisq(x, df, false, false);
    double cpp_p1_up_log = pchisq(x, df, false, true);

    double e_cpp_p1 = cdf_chisq(x, df, true);
    double e_cpp_p1_log = cdf_chisq_log(x, df, true);
    double e_cpp_p1_up = cdf_chisq(x, df, false);
    double e_cpp_p1_up_log = cdf_chisq_log(x, df, false);


    double eps = 1e-50;
    double eps_log = 1e-50;
    CHECK(cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));

    //now checking exposed, should be the same
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}
