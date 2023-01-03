#include<float.h>
#include<math.h>
#include<string.h>
#include<ostream>
#include<iomanip>
// for throwing error for wmw_test
#include <stdexcept>

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


/* log(sqrt(pi))== log(pi)/2 */
#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	
#endif


/* log(sqrt(2*pi)) == log(2*pi)/2 */
#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	
#endif


/* log(sqrt(pi/2))== log(pi/2)/2 */
#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	
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


#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif


//////////////////////////////////////////////////////////////////////////////
// END OF CONSTANTS
//////////////////////////////////////////////////////////////////////////////


//in: dpq.h
//interestingly, this macro MAY return a value, and may not 
//used for checking values are within bounds
#define R_P_bounds_01(x, x_min, x_max)	\
    if(x <= x_min) return R_DT_0;		\
    if(x >= x_max) return R_DT_1




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

/* C++ math header undefines any isnan macro. This file
   doesn't get C++ headers and so is safe. */
int R_isnancpp(double x)
{
    return (isnan(x) != 0);
}

static double myfmod(double x1, double x2)
{
    double q = x1 / x2;
    return x1 - floor(q) * x2;
}

double R_pow(double x, double y) /* = x ^ y */
{
    if(x == 1. || y == 0.)
	return(1.);
    if(x == 0.) {
	if(y > 0.) return(0.);
	/* y < 0 */return(ML_POSINF);
    }
    if (R_FINITE(x) && R_FINITE(y))
	return(pow(x,y));
    if (ISNAN(x) || ISNAN(y)) {
#ifdef IEEE_754
	return(x + y);
#else
	return(ML_NAN);
#endif
    }
    if(!R_FINITE(x)) {
	if(x > 0)		/* Inf ^ y */
	    return((y < 0.)? 0. : ML_POSINF);
	else {			/* (-Inf) ^ y */
	    if(R_FINITE(y) && y == floor(y)) /* (-Inf) ^ n */
		return((y < 0.) ? 0. : (myfmod(y,2.) ? x  : -x));
	}
    }
    if(!R_FINITE(y)) {
	if(x >= 0) {
	    if(y > 0)		/* y == +Inf */
		return((x >= 1)? ML_POSINF : 0.);
	    else		/* y == -Inf */
		return((x < 1) ? ML_POSINF : 0.);
	}
    }
    return(ML_NAN);		/* all other cases: (-Inf)^{+-Inf,
				   non-int}; (neg)^{+-Inf} */
}

double R_pow_di(double x, int n)
{
    double pow = 1.0;

    if (ISNAN(x)) return x;
    if (n != 0) {
	if (!R_FINITE(x)) return R_pow(x, (double)n);
	if (n < 0) { n = -n; x = 1/x; }
	for(;;) {
	    if(n & 01) pow *= x;
	    if(n >>= 1) x *= x; else break;
	}
    }
    return pow;
}

double NA_REAL = ML_NAN;
double R_PosInf = ML_POSINF, R_NegInf = ML_NEGINF;

//important helper function
#include <stdio.h>
#include <stdarg.h>
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




//from lgammacor.c
//Computes log(1+x)
double lgammacor(double x)
{
    const static double algmcs[15] = {
    +.1666389480451863247205729650822e+0,
    -.1384948176067563840732986059135e-4,
    +.9810825646924729426157171547487e-8,
    -.1809129475572494194263306266719e-10,
    +.6221098041892605227126015543416e-13,
    -.3399615005417721944303330599666e-15,
    +.2683181998482698748957538846666e-17,
    -.2868042435334643284144622399999e-19,
    +.3962837061046434803679306666666e-21,
    -.6831888753985766870111999999999e-23,
    +.1429227355942498147573333333333e-24,
    -.3547598158101070547199999999999e-26,
    +.1025680058010470912000000000000e-27,
    -.3401102254316748799999999999999e-29,
    +.1276642195630062933333333333333e-30
    };

    double tmp;

/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#define xmax  3.745194030963158e306

    if (x < 10)
    ML_ERR_return_NAN
    else if (x >= xmax) {
    ML_ERROR(ME_UNDERFLOW, "lgammacor");
    /* allow to underflow below */
    }
    else if (x < xbig) {
    tmp = 10 / x;
    return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
}



//sinpi
#ifdef HAVE_SINPI
#elif defined HAVE___SINPI
double sinpi(double x) {
    return __sinpi(x);
}
#else
// sin(pi * x)  -- exact when x = k/2  for all integer k
double sinpi(double x) {
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_ERR_return_NAN;

    x = fmod(x, 2.); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
    // map (-2,2) --> (-1,1] :
    if(x <= -1) x += 2.; else if (x > 1.) x -= 2.;
    if(x == 0. || x == 1.) return 0.;
    if(x ==  0.5)	return  1.;
    if(x == -0.5)	return -1.;
    // otherwise
    return sin(M_PI * x);
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


//gammalims, which is needed in gammafn
void gammalims(double *xmin, double *xmax2){
/* FIXME: Even better: If IEEE, #define these in nmath.h
	  and don't call gammalims() at all
*/
#ifdef IEEE_754
    *xmin = -170.5674972726612;
    *xmax2 =  171.61447887182298;/*(3 Intel/Sparc architectures)*/
#else
    double alnbig, alnsml, xln, xold;
    int i;

    alnsml = log(d1mach(1));
    *xmin = -alnsml;
    for (i=1; i<=10; ++i) {
	xold = *xmin;
	xln = log(*xmin);
	*xmin -= *xmin * ((*xmin + .5) * xln - *xmin - .2258 + alnsml) /
		(*xmin * xln + .5);
	if (fabs(*xmin - xold) < .005) {
	    *xmin = -(*xmin) + .01;
	    goto find_xmax2;
	}
    }

    /* unable to find xmin */

    ML_ERROR(ME_NOCONV, "gammalims");
    *xmin = *xmax2 = ML_NAN;

find_xmax2:

    alnbig = log(d1mach(2));
    *xmax2 = alnbig;
    for (i=1; i<=10; ++i) {
	xold = *xmax2;
	xln = log(*xmax2);
	*xmax2 -= *xmax2 * ((*xmax2 - .5) * xln - *xmax2 + .9189 - alnbig) /
		(*xmax2 * xln - .5);
	if (fabs(*xmax2 - xold) < .005) {
	    *xmax2 += -.01;
	    goto done;
	}
    }

    /* unable to find xmax2 */

    ML_ERROR(ME_NOCONV, "gammalims");
    *xmin = *xmax2 = ML_NAN;

done:
    *xmin = fmax2(*xmin, -(*xmax2) + 1);
#endif
}

//TODO: declaring lgammafn here, although only defined later
double lgammafn(double);


//stirlerr, which is needed in gammafn
/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 */

double stirlerr(double n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0) {
	nn = n + n;
	if (nn == (int)nn) return(sferr_halves[(int)nn]);
	return(lgammafn(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}




//gammanfn, which is needed in lgammafn_sign
double gammafn(double x)
{
    const static double gamcs[42] = {
	+.8571195590989331421920062399942e-2,
	+.4415381324841006757191315771652e-2,
	+.5685043681599363378632664588789e-1,
	-.4219835396418560501012500186624e-2,
	+.1326808181212460220584006796352e-2,
	-.1893024529798880432523947023886e-3,
	+.3606925327441245256578082217225e-4,
	-.6056761904460864218485548290365e-5,
	+.1055829546302283344731823509093e-5,
	-.1811967365542384048291855891166e-6,
	+.3117724964715322277790254593169e-7,
	-.5354219639019687140874081024347e-8,
	+.9193275519859588946887786825940e-9,
	-.1577941280288339761767423273953e-9,
	+.2707980622934954543266540433089e-10,
	-.4646818653825730144081661058933e-11,
	+.7973350192007419656460767175359e-12,
	-.1368078209830916025799499172309e-12,
	+.2347319486563800657233471771688e-13,
	-.4027432614949066932766570534699e-14,
	+.6910051747372100912138336975257e-15,
	-.1185584500221992907052387126192e-15,
	+.2034148542496373955201026051932e-16,
	-.3490054341717405849274012949108e-17,
	+.5987993856485305567135051066026e-18,
	-.1027378057872228074490069778431e-18,
	+.1762702816060529824942759660748e-19,
	-.3024320653735306260958772112042e-20,
	+.5188914660218397839717833550506e-21,
	-.8902770842456576692449251601066e-22,
	+.1527474068493342602274596891306e-22,
	-.2620731256187362900257328332799e-23,
	+.4496464047830538670331046570666e-24,
	-.7714712731336877911703901525333e-25,
	+.1323635453126044036486572714666e-25,
	-.2270999412942928816702313813333e-26,
	+.3896418998003991449320816639999e-27,
	-.6685198115125953327792127999999e-28,
	+.1146998663140024384347613866666e-28,
	-.1967938586345134677295103999999e-29,
	+.3376448816585338090334890666666e-30,
	-.5793070335782135784625493333333e-31
    };

    int i, n;
    double y;
    double sinpiy, value;

#ifdef NOMORE_FOR_THREADS
    static int ngam = 0;
    static double xmin = 0, xmax3 = 0., xsml = 0., dxrel2 = 0.;

    /* Initialize machine dependent constants, the first time gamma() is called.
	FIXME for threads ! */
    if (ngam == 0) {
	ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
	gammalims(&xmin, &xmax3);/*-> ./gammalims.c */
	xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
	/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
	dxrel2 = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 * (xmin, xmax3) are non-trivial, see ./gammalims.c
 * xsml = exp(.01)*DBL_MIN
 * dxrel2 = sqrt(DBL_EPSILON) = 2 ^ -26
*/
# define ngam 22
# define xmin -170.5674972726612
# define xmax3  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel2 1.490116119384765696e-8
#endif

    if(ISNAN(x)) return x;

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == round(x))) {
	ML_ERROR(ME_DOMAIN, "gammafn");
	return ML_NAN;
    }

    y = fabs(x);

    if (y <= 10) {

	/* Compute gamma(x) for -10 <= x <= 10
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */

	n = (int) x;
	if(x < 0) --n;
	y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
	--n;
	value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
	if (n == 0)
	    return value;/* x = 1.dddd = 1+y */

	if (n < 0) {
	    /* compute gamma(x) for -10 <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
	    if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel2) {
		ML_ERROR(ME_PRECISION, "gammafn");
	    }

	    /* The argument is so close to 0 that the result would overflow. */
	    if (y < xsml) {
		ML_ERROR(ME_RANGE, "gammafn");
		if(x > 0) return ML_POSINF;
		else return ML_NEGINF;
	    }

	    n = -n;

	    for (i = 0; i < n; i++) {
		value /= (x + i);
	    }
	    return value;
	}
	else {
	    /* gamma(x) for 2 <= x <= 10 */

	    for (i = 1; i <= n; i++) {
		value *= (y + i);
	    }
	    return value;
	}
    }
    else {
	/* gamma(x) for	 y = |x| > 10. */

	if (x > xmax3) {			/* Overflow */
	    ML_ERROR(ME_RANGE, "gammafn");
	    return ML_POSINF;
	}

	if (x < xmin) {			/* Underflow */
	    ML_ERROR(ME_UNDERFLOW, "gammafn");
	    return 0.;
	}

	if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
	    value = 1.;
	    for (i = 2; i < y; i++) value *= i;
	}
	else { /* normal case */
	    value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
			((2*y == (int)2*y)? stirlerr(y) : lgammacor(y)));
	}
	if (x > 0)
	    return value;

	if (fabs((x - (int)(x - 0.5))/x) < dxrel2){

	    /* The answer is less than half precision because */
	    /* the argument is too near a negative integer. */

	    ML_ERROR(ME_PRECISION, "gammafn");
	}

	sinpiy = sinpi(y);
	if (sinpiy == 0) {		/* Negative integer arg - overflow */
	    ML_ERROR(ME_RANGE, "gammafn");
	    return ML_POSINF;
	}

	return -M_PI / (y * sinpiy * value);
    }
}




//lgammafn
double lgammafn_sign(double x, int *sgn)
{
    double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static double xmax1 = 0.;
    static double dxrel = 0.;

    if (xmax1 == 0) {/* initialize machine dependent constants _ONCE_ */
	xmax1 = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
	dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax1  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
 */
#define xmax1  2.5327372760800758e+305
#define dxrel 1.490116119384765625e-8
#endif

    if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
    if(ISNAN(x)) return x;
#endif

    if (sgn != NULL && x < 0 && fmod(floor(-x), 2.) == 0)
	*sgn = -1;

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
	ML_ERROR(ME_RANGE, "lgamma");
	return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y < 1e-306) return -log(y); // denormalized range, R change
    if (y <= 10) return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax1) {
	ML_ERROR(ME_RANGE, "lgamma");
	return ML_POSINF;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
	if(x > 1e17)
	    return(x*(log(x) - 1.));
	else if(x > 4934720.)
	    return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
	else
#endif
	    return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sinpi(y));

    if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
	MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
	ML_ERR_return_NAN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

	/* The answer is less than half precision because
	 * the argument is too near a negative integer. */

	ML_ERROR(ME_PRECISION, "lgamma");
    }

    return ans;
}

double lgammafn(double x)
{
    return lgammafn_sign(x, NULL);
}


// From lbeta.c from:
double lbeta(double a, double b)
{
    double corr, p, q;

#ifdef IEEE_754
    if(ISNAN(a) || ISNAN(b))
    return a + b;
#endif
    p = q = a;
    if(b < p) p = b;/* := min(a,b) */
    if(b > q) q = b;/* := max(a,b) */

    /* both arguments must be >= 0 */
    if (p < 0)
    ML_ERR_return_NAN
    else if (p == 0) {
    return ML_POSINF;
    }
    else if (!R_FINITE(q)) { /* q == +Inf */
    return ML_NEGINF;
    }

    if (p >= 10) {
    /* p and q are big. */
    corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
    return log(q) * -0.5 + M_LN_SQRT_2PI + corr
    	+ (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
    }
    else if (q >= 10) {
    /* p is small, but q is big. */
    corr = lgammacor(q) - lgammacor(p + q);
    return lgammafn(p) + corr + p - p * log(p + q)
    	+ (q - 0.5) * log1p(-p / (p + q));
    }
    else {
    /* p and q are small: p <= q < 10. */
    /* R change for very small args */
    if (p < 1e-306) return lgamma(p) + (lgamma(q) - lgamma(p+q));
    else return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
    }
}




//bd0 from bd0.c
double bd0(double x, double np)
{
    double ej, s, s1, v;
    int j;

    if(!R_FINITE(x) || !R_FINITE(np) || np == 0.0) ML_ERR_return_NAN;

    if (fabs(x-np) < 0.1*(x+np)) {
	v = (x-np)/(x+np);  // might underflow to 0
	s = (x-np)*v;/* s using v -- change by MM */
	if(fabs(s) < DBL_MIN) return s;
	ej = 2*x*v;
	v = v*v;
	for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
					as |v| < .1,  v^2000 is "zero" */
	    ej *= v;// = v^(2j+1)
	    s1 = s+ej/((j<<1)+1);
	    if (s1 == s) /* last term was effectively 0 */
		return s1 ;
	    s = s1;
	}
    }
    /* else:  | x - np |  is not too small */
    return(x*log(x/np)+np-x);
}









//dnorm4 in dnorm.c
double dnorm4(double x, double mu, double sigma, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if(!R_FINITE(sigma)) return R_D__0;
    if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
    if (sigma <= 0) {
	if (sigma < 0) ML_ERR_return_NAN;
	/* sigma == 0 */
	return (x == mu) ? ML_POSINF : R_D__0;
    }
    x = (x - mu) / sigma;

    if(!R_FINITE(x)) return R_D__0;

    x = fabs (x);
    if (x >= 2 * sqrt(DBL_MAX)) return R_D__0;
    if (give_log)
	return -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma));
    //  M_1_SQRT_2PI = 1 / sqrt(2 * pi)
#ifdef MATHLIB_FAST_dnorm
    // and for R <= 3.0.x and R-devel upto 2014-01-01:
    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
#else
    // more accurate, less fast :
    if (x < 5)    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;

    /* ELSE:

     * x*x  may lose upto about two digits accuracy for "large" x
     * Morten Welinder's proposal for PR#15620
     * https://bugs.r-project.org/bugzilla/show_bug.cgi?id=15620

     * -- 1 --  No hoop jumping when we underflow to zero anyway:

     *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
     *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
     * but "thanks" to denormalized numbers, underflow happens a bit later,
     *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
     * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
     * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
     *              =IEEE=  38.58601
     * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
     */
    if (x > sqrt(-2*M_LN2*(DBL_MIN_EXP + 1-DBL_MANT_DIG))) return 0.;

    /* Now, to get full accurary, split x into two parts,
     *  x = x1+x2, such that |x2| <= 2^-16.
     * Assuming that we are using IEEE doubles, that means that
     * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).

     * If we do not have IEEE this is still an improvement over the naive formula.
     */
    double x1 = //  R_forceint(x * 65536) / 65536 =
	ldexp( R_forceint(ldexp(x, 16)), -16);
    double x2 = x - x1;
    return M_1_SQRT_2PI / sigma *
	(exp(-0.5 * x1 * x1) * exp( (-0.5 * x2 - x1) * x2 ) );
#endif
}


//synonym for dnorm
//redefinition
#define dnorm dnorm4



//dpnorm from pgamma.c
/*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *		 dnorm (x, 0, 1, FALSE)
 *	   ----------------------------------
 *	   pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
static double
dpnorm (double x, int lower_tail, double lp){
    /*
     * So as not to repeat a pnorm call, we expect
     *
     *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
     *
     * but use it only in the non-critical case where either x is small
     * or p==exp(lp) is close to 1.
     */

    if (x < 0) {
	x = -x;
	lower_tail = !lower_tail;
    }

    if (x > 10 && !lower_tail) {
	double term = 1 / x;
	double sum = term;
	double x2 = x * x;
	double i = 1;

	do {
	    term *= -i / x2;
	    sum += term;
	    i += 2;
	} while (fabs (term) > DBL_EPSILON * sum);

	return 1 / sum;
    } else {
	//FIXME double d = dnorm (x, 0., 1., FALSE);
    //using Rboolean would be fine, but not here...
	double d = dnorm (x, 0., 1., false);
	return d / exp (lp);
    }
}


//logcf from pgamma.c
/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxilary in log1pmx() and lgamma1p()
 */
static double
logcf (double x, double i, double d,
       double eps /* ~ relative tolerance */)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
	double c3 = c2*c2*x;
	c2 += d;
	c4 += d;
	a1 = c4 * a2 - c3 * a1;
	b1 = c4 * b2 - c3 * b1;

	c3 = c1 * c1 * x;
	c1 += d;
	c4 += d;
	a2 = c4 * a1 - c3 * a2;
	b2 = c4 * b1 - c3 * b2;

	if (fabs (b2) > scalefactor) {
	    a1 /= scalefactor;
	    b1 /= scalefactor;
	    a2 /= scalefactor;
	    b2 /= scalefactor;
	} else if (fabs (b2) < 1 / scalefactor) {
	    a1 *= scalefactor;
	    b1 *= scalefactor;
	    a2 *= scalefactor;
	    b2 *= scalefactor;
	}
    }

    return a2 / b2;
}





//log1pmx from pgamma.c
/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double log1pmx (double x)
{
    static const double minLog1Value = -0.79149064;

    if (x > 1 || x < minLog1Value)
	return log1p(x) - x;
    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
	    * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
	    * ---------------------------------------------
	    * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
	   */
	double r = x / (2 + x), y = r * r;
	if (fabs(x) < 1e-2) {
	    static const double two = 2;
	    return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
			    two / 3) * y - x);
	} else {
	    static const double tol_logcf = 1e-14;
	    return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
	}
    }
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



//lgamma1p
/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p (double a)
{
    const double eulers_const =	 0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const int N = 40;
    static const double coeffs[40] = {
	0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
	0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
	0.2058080842778454787900092413529198e-1,
	0.7385551028673985266273097291406834e-2,
	0.2890510330741523285752988298486755e-2,
	0.1192753911703260977113935692828109e-2,
	0.5096695247430424223356548135815582e-3,
	0.2231547584535793797614188036013401e-3,
	0.9945751278180853371459589003190170e-4,
	0.4492623673813314170020750240635786e-4,
	0.2050721277567069155316650397830591e-4,
	0.9439488275268395903987425104415055e-5,
	0.4374866789907487804181793223952411e-5,
	0.2039215753801366236781900709670839e-5,
	0.9551412130407419832857179772951265e-6,
	0.4492469198764566043294290331193655e-6,
	0.2120718480555466586923135901077628e-6,
	0.1004322482396809960872083050053344e-6,
	0.4769810169363980565760193417246730e-7,
	0.2271109460894316491031998116062124e-7,
	0.1083865921489695409107491757968159e-7,
	0.5183475041970046655121248647057669e-8,
	0.2483674543802478317185008663991718e-8,
	0.1192140140586091207442548202774640e-8,
	0.5731367241678862013330194857961011e-9,
	0.2759522885124233145178149692816341e-9,
	0.1330476437424448948149715720858008e-9,
	0.6422964563838100022082448087644648e-10,
	0.3104424774732227276239215783404066e-10,
	0.1502138408075414217093301048780668e-10,
	0.7275974480239079662504549924814047e-11,
	0.3527742476575915083615072228655483e-11,
	0.1711991790559617908601084114443031e-11,
	0.8315385841420284819798357793954418e-12,
	0.4042200525289440065536008957032895e-12,
	0.1966475631096616490411045679010286e-12,
	0.9573630387838555763782200936508615e-13,
	0.4664076026428374224576492565974577e-13,
	0.2273736960065972320633279596737272e-13,
	0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };

    const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    const double tol_logcf = 1e-14;
    double lgam;
    int i;

    if (fabs (a) >= 0.5)
	return lgammafn (a + 1);

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
    for (i = N - 1; i >= 0; i--)
	lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx (a);
} /* lgamma1p */




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



//////////////////////////////////////////////////////////////////////////////
//Now for pgamma functions
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//Elementary functions building up to pt
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
//gam elementary functions from toms708.c
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
//More elementary functions
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
// Inverse t cdf - quantile function
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Inverse gamma cdf
//////////////////////////////////////////////////////////////////////////////

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
// q t c o d e
//////////////////////////////////////////////////////////////////////////////


#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */

double tanpi(double x) {
    return __tanpi(x);
}





double fmin2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? x : y;
}








//////////////////////////////////////////////////////////////////////////////
// Wilcox
//////////////////////////////////////////////////////////////////////////////

#define WILCOX_MAX 50

int imax2(int x, int y)
{
    return (x < y) ? y : x;
}

static double ***w; /* to store  cwilcox(i,j,k) -> w[i][j][k] */
static int allocated_m, allocated_n;

static void
w_free(int m, int n)
{
    int i, j;

    for (i = m; i >= 0; i--) {
	for (j = n; j >= 0; j--) {
	    if (w[i][j] != 0)
		free((void *) w[i][j]);
	}
	free((void *) w[i]);
    }
    free((void *) w);
    w = 0; allocated_m = allocated_n = 0;
}

static void
w_init_maybe(int m, int n)
{
    int i;

    if (m > n) {
	i = n; n = m; m = i;
    }
    if (w && (m > allocated_m || n > allocated_n))
	w_free(allocated_m, allocated_n); /* zeroes w */

    if (!w) { /* initialize w[][] */
	m = imax2(m, WILCOX_MAX);
	n = imax2(n, WILCOX_MAX);
	w = (double ***) calloc((size_t) m + 1, sizeof(double **));
#ifdef MATHLIB_STANDALONE
	if (!w) MATHLIB_ERROR(_("wilcox allocation error %d"), 1);
#endif
	for (i = 0; i <= m; i++) {
	    w[i] = (double **) calloc((size_t) n + 1, sizeof(double *));
#ifdef MATHLIB_STANDALONE
	    /* the apparent leak here in the in-R case should be
	       swept up by the on.exit action */
	    if (!w[i]) {
		/* first free all earlier allocations */
		w_free(i-1, n);
		MATHLIB_ERROR(_("wilcox allocation error %d"), 2);
	    }
#endif
	}
	allocated_m = m; allocated_n = n;
    }
}

static void
w_free_maybe(int m, int n)
{
    if (m > WILCOX_MAX || n > WILCOX_MAX)
	w_free(m, n);
}


# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
// # define MATHLIB_WARNING2(fmt,x,x2)	printf(fmt,x,x2)

# define R_nonint(x) 	  (fabs((x) - R_forceint(x)) > 1e-7*fmax2(1., fabs(x)))

/* This counts the number of choices with statistic = k */
static double
cwilcox(int k, int m, int n)
{
    int c, u, i, j, l;

// NOTE XXX Removing this, because it is a pain to need to define all the errors
// #ifndef MATHLIB_STANDALONE
//     R_CheckUserInterrupt();
//
//     MATHLIB_WARNING("Warning: R_CheckUserInterrupt in cwilcox (%g)\n", 0)
// #endif

    u = m * n;
    if (k < 0 || k > u)
	return(0);
    c = (int)(u / 2);
    if (k > c)
	k = u - k; /* hence  k <= floor(u / 2) */
    if (m < n) {
	i = m; j = n;
    } else {
	i = n; j = m;
    } /* hence  i <= j */

    if (j == 0) /* and hence i == 0 */
	return (k == 0);


    /* We can simplify things if k is small.  Consider the Mann-Whitney 
       definition, and sort y.  Then if the statistic is k, no more 
       than k of the y's can be <= any x[i], and since they are sorted 
       these can only be in the first k.  So the count is the same as
       if there were just k y's. 
    */
    if (j > 0 && k < j) return cwilcox(k, i, k);    
    
    if (w[i][j] == 0) {
	w[i][j] = (double *) calloc((size_t) c + 1, sizeof(double));
#ifdef MATHLIB_STANDALONE
	if (!w[i][j]) MATHLIB_ERROR(_("wilcox allocation error %d"), 3);
#endif
	for (l = 0; l <= c; l++)
	    w[i][j][l] = -1;
    }
    if (w[i][j][k] < 0) {
	if (j == 0) /* and hence i == 0 */
	    w[i][j][k] = (k == 0);
	else
	    w[i][j][k] = cwilcox(k - j, i - 1, j) + cwilcox(k, i, j - 1);

    }
    return(w[i][j][k]);
}

//choose


#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_visible __attribute__ ((visibility ("default")))
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_visible
# define attribute_hidden
#endif


double attribute_hidden lfastchoose(double n, double k)
{
    return -log(n + 1.) - lbeta(n - k + 1., k + 1.);
}
/* mathematically the same:
   less stable typically, but useful if n-k+1 < 0 : */
static
double lfastchoose2(double n, double k, int *s_choose)
{
    double r;
    r = lgammafn_sign(n - k + 1., s_choose);
    return lgammafn(n + 1.) - lgammafn(k + 1.) - r;
}

#define ODD(_K_) ((_K_) != 2 * floor((_K_) / 2.))

#define R_IS_INT(x)  (!R_nonint(x))

double lchoose(double n, double k)
{
    double k0 = k;
    k = R_forceint(k);
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if(ISNAN(n) || ISNAN(k)) return n + k;
#endif
//NOTE XXX removing this
// #ifndef MATHLIB_STANDALONE
//     R_CheckStack();
// #endif
    if (fabs(k - k0) > 1e-7)
	MATHLIB_WARNING2(_("'k' (%.2f) must be integer, rounded to %.0f"), k0, k);
    if (k < 2) {
	if (k <	 0) return ML_NEGINF;
	if (k == 0) return 0.;
	/* else: k == 1 */
	return log(fabs(n));
    }
    /* else: k >= 2 */
    if (n < 0) {
	return lchoose(-n+ k-1, k);
    }
    else if (R_IS_INT(n)) {
	n = R_forceint(n);
	if(n < k) return ML_NEGINF;
	/* k <= n :*/
	if(n - k < 2) return lchoose(n, n-k); /* <- Symmetry */
	/* else: n >= k+2 */
	return lfastchoose(n, k);
    }
    /* else non-integer n >= 0 : */
    if (n < k-1) {
	int s;
	return lfastchoose2(n, k, &s);
    }
    return lfastchoose(n, k);
}

#define k_small_max 30
/* 30 is somewhat arbitrary: it is on the *safe* side:
 * both speed and precision are clearly improved for k < 30.
*/
double choose(double n, double k)
{
    double r, k0 = k;
    k = R_forceint(k);
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if(ISNAN(n) || ISNAN(k)) return n + k;
#endif

//NOTE XXX avoiding MATHLIB_STANADLONE to avoid checking stack and defining SEXP
// #ifndef MATHLIB_STANDALONE
//     R_CheckStack();
// #endif
    if (fabs(k - k0) > 1e-7)
	MATHLIB_WARNING2(_("'k' (%.2f) must be integer, rounded to %.0f"), k0, k);
    if (k < k_small_max) {
	int j;
	if(n-k < k && n >= 0 && R_IS_INT(n))
	    k = R_forceint(n-k); /* <- Symmetry, ensure k still integer */
	if (k <	 0) return 0.;
	if (k == 0) return 1.;
	/* else: k >= 1 */
	r = n;
	for(j = 2; j <= k; j++)
	    r *= (n-j+1)/j;
	return R_IS_INT(n) ? R_forceint(r) : r;
	/* might have got rounding errors */
    }
    /* else: k >= k_small_max */
    if (n < 0) {
	r = choose(-n+ k-1, k);
	if (ODD(k)) r = -r;
	return r;
    }
    else if (R_IS_INT(n)) {
	n = R_forceint(n);
	if(n < k) return 0.;
	if(n - k < k_small_max) return choose(n, n-k); /* <- Symmetry */
	return R_forceint(exp(lfastchoose(n, k)));
    }
    /* else non-integer n >= 0 : */
    if (n < k-1) {
	int s_choose;
	r = lfastchoose2(n, k, /* -> */ &s_choose);
	return s_choose * exp(r);
    }
    return exp(lfastchoose(n, k));
}



// pwilcox
//


#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */

#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */

#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))

/* args have the same meaning as R function pwilcox */
double pwilcox(double q, double m, double n, int lower_tail, int log_p)
{
    int i;
    double c, p;

#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(m) || ISNAN(n))
	return(q + m + n);
#endif
    if (!R_FINITE(m) || !R_FINITE(n))
	ML_WARN_return_NAN;
    m = R_forceint(m);
    n = R_forceint(n);
    if (m <= 0 || n <= 0)
	ML_WARN_return_NAN;

    q = floor(q + 1e-7);

    if (q < 0.0)
	return(R_DT_0);
    if (q >= m * n)
	return(R_DT_1);

    int mm = (int) m, nn = (int) n;
    w_init_maybe(mm, nn);
    c = choose(m + n, n);
    p = 0;
    /* Use summation of probs over the shorter range */
    if (q <= (m * n / 2)) {
	for (i = 0; i <= q; i++)
	    p += cwilcox(i, mm, nn) / c;
    }
    else {
	q = m * n - q;
	for (i = 0; i < q; i++)
	    p += cwilcox(i, mm, nn) / c;
	lower_tail = !lower_tail; /* p = 1 - p; */
    }

    return(R_DT_val(p));
} /* pwilcox */


//////////////////////////////////////////////////////////////////////////////
// Sign function
//////////////////////////////////////////////////////////////////////////////

/**
 * From sign.c
 *    This function computes the  'signum(.)' function:
 *
 *    	sign(x) =  1  if x > 0
 *    	sign(x) =  0  if x == 0
 *    	sign(x) = -1  if x < 0
 */

double sign(double x)
{
    if (ISNAN(x))
	return x;
    return ((x > 0) ? 1 : ((x == 0)? 0 : -1));
}



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


//cdf of wilcoxon distribution, log_p=false
double cdf_wilcoxon(double q, double m, double n, int lower_tail){
    return pwilcox(q, m, n, lower_tail, false);
}

//cdf of wilcoxon distribution, log_p=true
double cdf_wilcoxon_log(double q, double m, double n, int lower_tail){
    return pwilcox(q, m, n, lower_tail, true);
}



//////////////////////////////////////////////////////////////////////////////
