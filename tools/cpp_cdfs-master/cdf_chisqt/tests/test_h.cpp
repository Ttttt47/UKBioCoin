#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../../doctest/doctest.h"

//include others
#include <iostream>
#include "../cdf_base.h"


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
    double r_p1_up_log = -1.279427536176580559867943520657718181610107421875000000000000;


    double e_cpp_p1 = cdf_t(x, df, true);
    double e_cpp_p1_log = cdf_t_log(x, df, true);
    double e_cpp_p1_up = cdf_t(x, df, false);
    double e_cpp_p1_up_log = cdf_t_log(x, df, false);
//     int prec = 100;
//     std::cout << std::setprecision(prec) << r_p1 << std::endl;
//     std::cout << std::setprecision(prec) << cpp_p1 << std::endl;
    
    double eps = 1e-50;
    double eps_log = 1e-50;

    //now checking exposed
    CHECK(e_cpp_p1 == doctest::Approx(r_p1).epsilon(eps));
    CHECK(e_cpp_p1_log == doctest::Approx(r_p1_log).epsilon(eps_log));
    CHECK(e_cpp_p1_up == doctest::Approx(r_p1_up).epsilon(eps));
    CHECK(e_cpp_p1_up_log == doctest::Approx(r_p1_up_log).epsilon(eps_log));
}//testing pt



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
    double r_p1_up_log = -0.270000000000000017763568394002504646778106689453125000000000;

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
