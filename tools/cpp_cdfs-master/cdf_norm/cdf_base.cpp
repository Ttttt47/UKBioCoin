#include<float.h>
#include<math.h>
#include<string.h>
#include<ostream>
#include<iomanip>
// for throwing error for wmw_test
#include <stdexcept>
#include <stdio.h>
#include <stdarg.h>

/*
 * All of this code has been taken from the R Core library with
 * minor modifications.
 *
 * In this file, allmost all the comments have been removed.
 *  See cdf_base_with_comments.cpp for all the comments.
 *
 *
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2001 The R Core Team
 *  and
 *	Catherine Loader, catherine@research.bell-labs.com.
 *	October 23, 2000.
 *
 *  bratio written by Alfred H. Morris, Jr.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

//from src/extra/tre/regerror.c
#define gettext(s) s

#define _(STRING) gettext(STRING)

#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	((-1.0) / 0.0)
#define ML_NAN		(0.0 / 0.0)

//in: dpq.h
#define give_log log_p
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

// #define FALSE false
// #define TRUE true

//from nmath.h
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ME_RANGE	2
#define ME_DOMAIN	1
#define ME_NOCONV	4

//in pgamma.c
/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

//isnan is cpp function
#define ISNAN(x) (isnan(x)!=0)

//min and max from toms708.c
#undef min
#define min(a,b) ((a < b)?a:b)
#undef max
#define max(a,b) ((a > b)?a:b)


#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif


/* log(sqrt(2*pi)) == log(2*pi)/2 */
#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	
#endif


#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif


#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif


#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

//////////////////////////////////////////////////////////////////////////////
// END OF CONSTANTS
//////////////////////////////////////////////////////////////////////////////



// Arith.h defines it
//from nmath,h
#ifndef R_FINITE
#ifdef HAVE_WORKING_ISFINITE
/* isfinite is defined in <math.h> according to C99 */
# define R_FINITE(x)    isfinite(x)
#else
# define R_FINITE(x)    R_finite(x)
#endif
#endif

//from mlutils.h
int R_finite(double x)
{
#ifdef HAVE_WORKING_ISFINITE
    return isfinite(x);
# else
    return (!isnan(x) & (x != ML_POSINF) & (x != ML_NEGINF));
#endif
}



//from errors.c
//TODO: fix this warning issue in order for chebyshev_eval
// and other functions to work properly and be able to throw errors
void warning(const std::string format, ...)
{
//     char buf[BUFSIZE], *p;
//     RCNTXT *c = R_GlobalContext;
// 
//     va_list(ap);
//     va_start(ap, format);
//     Rvsnprintf(buf, min(BUFSIZE, R_WarnLength+1), format, ap);
//     va_end(ap);
//     p = buf + strlen(buf) - 1;
//     if(strlen(buf) > 0 && *p == '\n') *p = '\0';
//     RprintTrunc(buf);
//     if (c && (c->callflag & CTXT_BUILTIN)) c = c->nextcontext;
//     if (c == R_GlobalContext && R_BCIntActive)
//         warningcall(R_getBCInterpreterExpression(), "%s", buf);
//     else
//         warningcall(c ? c->call : R_NilValue, "%s", buf);

    //do nothing
    std::cerr << "warning:" << format;
}

//in nmath.h
//
// #define MATHLIB_ERROR(fmt,x)		error(fmt,x);
#define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
#define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
#define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
#define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)
#define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) warning(fmt,x,x2,x3,x4,x5)
 

#define ML_ERROR(x, s) { \
   if(x > ME_DOMAIN) { \
       std::string msg = ""; \
       switch(x) { \
       case ME_DOMAIN: \
	   msg = ("argument out of domain in '%s'\n");	\
	   break; \
       case ME_RANGE: \
	   msg = ("value out of range in '%s'\n");	\
	   break; \
       case ME_NOCONV: \
	   msg = ("convergence failed in '%s'\n");	\
	   break; \
       case ME_PRECISION: \
	   msg = ("full precision may not have been achieved in '%s'\n"); \
	   break; \
       case ME_UNDERFLOW: \
	   msg = ("underflow occurred in '%s'\n");	\
	   break; \
       } \
       MATHLIB_WARNING(msg, s); \
   } \
}

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }

/* moved from dpq.h */
#ifdef HAVE_NEARYINT
# define R_forceint(x)   nearbyint()
#else
# define R_forceint(x)   round(x)
#endif

//in dpq.h
#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */

#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))

#ifdef DEBUG_bratio
# define R_ifDEBUG_printf(...) REprintf(__VA_ARGS__)
#else
# define R_ifDEBUG_printf(...)
#endif


//helps to define REprintf
#ifdef MATHLIB_STANDALONE


double NA_REAL = ML_NAN;
double R_PosInf = ML_POSINF, R_NegInf = ML_NEGINF;

//important helper function
void attribute_hidden REprintf(const char *format, ...)
{
    va_list(ap);
    va_start(ap, format);
    fprintf(stderr, format, ap);
    va_end(ap);
}

#endif /* MATHLIB_STANDALONE */



//Rboolean defined in 
//src/include/R_ext/Boolean.h
//annoyingly, the conversion from boolean to enum will not work easily...
//so just renaming Rboolean as alias for bool
#ifndef R_EXT_BOOLEAN_H_
#define R_EXT_BOOLEAN_H_

#undef FALSE
#undef TRUE

#ifdef  __cplusplus
extern "C" {
#endif
// typedef enum { FALSE = 0, TRUE /*, MAYBE */ } Rboolean;
typedef bool Rboolean;

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_BOOLEAN_H_ */


//////////////////////////////////////////////////////////////////////////////
// ELEMENTARY FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

// attribute_hidden double Rf_d1mach(int i)
double Rf_d1mach(int i){
    switch(i) {
    case 1: return DBL_MIN;
    case 2: return DBL_MAX;

    case 3: /* = FLT_RADIX  ^ - DBL_MANT_DIG
	      for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON */
	return 0.5*DBL_EPSILON;

    case 4: /* = FLT_RADIX  ^ (1- DBL_MANT_DIG) =
	      for IEEE:  = 2^-52 = DBL_EPSILON */
	return DBL_EPSILON;

    case 5: return M_LOG10_2;

    default: return 0.0;
    }
}


double d1mach(int *i){
    return Rf_d1mach(*i);
}

//overloading
double d1mach(int i){
    return Rf_d1mach(i);
}


double chebyshev_eval(double x, const double *a, const int n)
{
    double b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) ML_ERR_return_NAN;

    if (x < -1.1 || x > 1.1) ML_ERR_return_NAN;

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

int chebyshev_init(double *dos, int nos, double eta){
    int i, ii;
    double err;

    if (nos < 1)
    return 0;

    err = 0.0;
    i = 0;			/* just to avoid compiler warnings */
    for (ii=1; ii<=nos; ii++) {
    i = nos - ii;
    err += fabs(dos[i]);
    if (err > eta) {
        return i;
    }
    }
    return i;
}


//from log1p.c
//Computes log(1+x)
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
// #include "nmath.h"

/* want to compile log1p as Rlog1p if HAVE_LOG1P && !HAVE_WORKING_LOG1P */
#if defined(HAVE_LOG1P) && !defined(HAVE_WORKING_LOG1P)
#undef HAVE_LOG1P
#endif

#ifndef HAVE_LOG1P
double log1p(double x)
{
    /* series for log1p on the interval -.375 to .375
     *				     with weighted error   6.35e-32
     *				      log weighted error  31.20
     *			    significant figures required  30.93
     *				 decimal places required  32.01
     */
    const static double alnrcs[43] = {
    +.10378693562743769800686267719098e+1,
    -.13364301504908918098766041553133e+0,
    +.19408249135520563357926199374750e-1,
    -.30107551127535777690376537776592e-2,
    +.48694614797154850090456366509137e-3,
    -.81054881893175356066809943008622e-4,
    +.13778847799559524782938251496059e-4,
    -.23802210894358970251369992914935e-5,
    +.41640416213865183476391859901989e-6,
    -.73595828378075994984266837031998e-7,
    +.13117611876241674949152294345011e-7,
    -.23546709317742425136696092330175e-8,
    +.42522773276034997775638052962567e-9,
    -.77190894134840796826108107493300e-10,
    +.14075746481359069909215356472191e-10,
    -.25769072058024680627537078627584e-11,
    +.47342406666294421849154395005938e-12,
    -.87249012674742641745301263292675e-13,
    +.16124614902740551465739833119115e-13,
    -.29875652015665773006710792416815e-14,
    +.55480701209082887983041321697279e-15,
    -.10324619158271569595141333961932e-15,
    +.19250239203049851177878503244868e-16,
    -.35955073465265150011189707844266e-17,
    +.67264542537876857892194574226773e-18,
    -.12602624168735219252082425637546e-18,
    +.23644884408606210044916158955519e-19,
    -.44419377050807936898878389179733e-20,
    +.83546594464034259016241293994666e-21,
    -.15731559416479562574899253521066e-21,
    +.29653128740247422686154369706666e-22,
    -.55949583481815947292156013226666e-23,
    +.10566354268835681048187284138666e-23,
    -.19972483680670204548314999466666e-24,
    +.37782977818839361421049855999999e-25,
    -.71531586889081740345038165333333e-26,
    +.13552488463674213646502024533333e-26,
    -.25694673048487567430079829333333e-27,
    +.48747756066216949076459519999999e-28,
    -.92542112530849715321132373333333e-29,
    +.17578597841760239233269760000000e-29,
    -.33410026677731010351377066666666e-30,
    +.63533936180236187354180266666666e-31,
    };

#ifdef NOMORE_FOR_THREADS
    static int nlnrel = 0;
    static double xmin = 0.0;

    if (xmin == 0.0) xmin = -1 + sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)); */
    if (nlnrel == 0) /* initialize chebychev coefficients */
    nlnrel = chebyshev_init(alnrcs, 43, DBL_EPSILON/20);/*was .1*d1mach(3)*/
#else
# define nlnrel 22
    const static double xmin = -0.999999985;
/* 22: for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
#endif

    if (x == 0.) return 0.;/* speed */
    if (x == -1) return(ML_NEGINF);
    if (x  < -1) ML_ERR_return_NAN;

    if (fabs(x) <= .375) {
        /* Improve on speed (only);
       again give result accurate to IEEE double precision: */
    if(fabs(x) < .5 * DBL_EPSILON)
        return x;

    if( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
        return x * (1 - .5 * x);
    /* else */
    return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
    }
    /* else */
    if (x < xmin) {
    /* answer less than half precision because x too near -1 */
    ML_ERROR(ME_PRECISION, "log1p");
    }
    return log(1 + x);
}
#endif



//fmax2
double fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}


//expm1 from expm1.c
double expm1(double x){
    double y, a = fabs(x);

    if (a < DBL_EPSILON) return x;
    if (a > 0.697) return exp(x) - 1;  /* negligible cancellation */

    if (a > 1e-8)
	y = exp(x) - 1;
    else /* Taylor expansion, more accurate in this range */
	y = (x / 2 + 1) * x;

    /* Newton step for solving   log(1 + y) = x   for y : */
    /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
    y -= (1 + y) * (log1p (y) - x);
    return y;
}




//////////////////////////////////////////////////////////////////////////////
//Another round of functions, building to ppois_asymp
//////////////////////////////////////////////////////////////////////////////


#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))


//pnorm_both from pnorm.c
#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p){
/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/
    const static double a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    };
    const static double b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
    };
    const static double c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
    };
    const static double d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
    };
    const static double p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
    };
    const static double q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
    };

    double xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
    double min = DBL_MIN;
#endif
    int i, lower, upper;

#ifdef IEEE_754
    if(ISNAN(x)) { *cum = *ccum = x; return; }
#endif

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	    xsq = x * x;
	    xnum = a[4] * xsq;
	    xden = xsq;
	    for (i = 0; i < 3; ++i) {
		xnum = (xnum + a[i]) * xsq;
		xden = (xden + b[i]) * xsq;
	    }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  *cum = 0.5 + temp;
	if(upper) *ccum = 0.5 - temp;
	if(log_p) {
	    if(lower)  *cum = log(*cum);
	    if(upper) *ccum = log(*ccum);
	}
    }
    else if (y <= M_SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	    xnum = (xnum + c[i]) * y;
	    xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = trunc(X * SIXTEN) / SIXTEN;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	    *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
	    if((lower && x > 0.) || (upper && x <= 0.))			\
		  *ccum = log1p(-exp(-xsq * xsq * 0.5) *		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    *ccum = 1.0 - *cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	}

	do_del(y);
	swap_tail;
    }

/* else	  |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !

 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((log_p && y < 1e170) /* avoid underflow below */
	/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

	 xsq = x*x;

	 if(xsq * DBL_EPSILON < 1.)
	    del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	 else
	    del = 0.;
	 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
	 *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

 	 swap_tail;

	 [Yes, but xsq might be infinite.]

	*/
	    || (lower && -37.5193 < x  &&  x < 8.2924)
	    || (upper && -8.2924  < x  &&  x < 37.5193)
	) {

	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	    xnum = (xnum + p[i]) * xsq;
	    xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (M_1_SQRT_2PI - temp) / y;

	do_del(x);
	swap_tail;
    } else { /* large x such that probs are 0 or 1 */
	if(x > 0) {	*cum = R_D__1; *ccum = R_D__0;	}
	else {	        *cum = R_D__0; *ccum = R_D__1;	}
    }


#ifdef NO_DENORMS
    /* do not return "denormalized" -- we do in R */
    if(log_p) {
	if(*cum > -min)	 *cum = -0.;
	if(*ccum > -min)*ccum = -0.;
    }
    else {
	if(*cum < min)	 *cum = 0.;
	if(*ccum < min)	*ccum = 0.;
    }
#endif
    return;
}






//pnorm5 from pnorm.c
double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p){
    double p, cp;

    /* Note: The structure of these checks has been carefully thought through.
     * For example, if x == mu and sigma == 0, we get the correct answer 1.
     */
#ifdef IEEE_754
    if(ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
    if (sigma <= 0) {
	if(sigma < 0) ML_ERR_return_NAN;
	/* sigma = 0 : */
	return (x < mu) ? R_DT_0 : R_DT_1;
    }
    p = (x - mu) / sigma;
    if(!R_FINITE(p))
	return (x < mu) ? R_DT_0 : R_DT_1;
    x = p;

    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

    return(lower_tail ? p : cp);
}


//need the synonym
#define pnorm pnorm5




//NOTE: minor change here from char* to string
//and changed _("X") to ("X"), i.e. no underscore

#define ML_WARNING(x, s) { \
   if(x > ME_DOMAIN) { \
       std::string msg = ""; \
       switch(x) { \
       case ME_DOMAIN: \
	   msg = ("argument out of domain in '%s'\n");	\
	   break; \
       case ME_RANGE: \
	   msg = ("value out of range in '%s'\n");	\
	   break; \
       case ME_NOCONV: \
	   msg = ("convergence failed in '%s'\n");	\
	   break; \
       case ME_PRECISION: \
	   msg = ("full precision may not have been achieved in '%s'\n"); \
	   break; \
       case ME_UNDERFLOW: \
	   msg = ("underflow occurred in '%s'\n");	\
	   break; \
       } \
       MATHLIB_WARNING(msg, s); \
   } \
}


#define ML_WARN_return_NAN { ML_WARNING(ME_DOMAIN, ""); return ML_NAN; }

#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_WARN_return_NAN


#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */

#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */


#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    ML_WARN_return_NAN;				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    ML_WARN_return_NAN;				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }


#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */

#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */


#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/

// #define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */

#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))




double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p)
{
    double p_, q, r, val;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
	return p + mu + sigma;
#endif
    R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

    if(sigma  < 0)	ML_WARN_return_NAN;
    if(sigma == 0)	return mu;

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

#ifdef DEBUG_qnorm
    REprintf("qnorm(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
	     p,mu,sigma, lower_tail, log_p, q);
#endif


/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* |p~ - 0.5| <= .425  <==> 0.075 <= p~ <= 0.925 */
        r = .180625 - q * q; // = .425^2 - q^2  >= 0
	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary :
	    *  r := log(p~);  p~ = min(p, 1-p) < 0.075 :  */
	if(log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0))) {
	    r = p;
	} else {
	    r = log( (q > 0) ? R_DT_CIv(p) /* 1-p */ : p_ /* = R_DT_Iv(p) ^=  p */);
	}
	// r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 ) :
        r = sqrt(-r);
#ifdef DEBUG_qnorm
	REprintf("\t close to 0 or 1: r = %7g\n", r);
#endif
        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else if(r >= 816) { // p is *extremly* close to 0 or 1 - only possibly when log_p =TRUE
	    // Using the asymptotical formula -- is *not* optimal but uniformly better than branch below
	    val = r * M_SQRT2;
        }
	else { // p is very close to  0 or 1:  r > 5 <==> min(p,1-p) < exp(-25) = 1.3888..e-11
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

	if(q < 0.0)
	    val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

#define qnorm qnorm5





//////////////////////////////////////////////////////////////////////////////
// exposed functions
//////////////////////////////////////////////////////////////////////////////



//cdf of normal dist, no log
double cdf_norm(double x, double mu, double sigma, int lower_tail){
    return pnorm(x, mu, sigma, lower_tail, false);
}

//cdf of normal dist, log
double cdf_norm_log(double x, double mu, double sigma, int lower_tail){
    return pnorm(x, mu, sigma, lower_tail, true);
}

//quantile of normal dist, no log
double quantile_norm(double p, double mu, double sigma, int lower_tail){
    return qnorm(p, mu, sigma, lower_tail, false);
}

//quantile of normal dist, log
double quantile_norm_log(double p, double mu, double sigma, int lower_tail){
    return qnorm(p, mu, sigma, lower_tail, true);
}

//////////////////////////////////////////////////////////////////////////////
