#include<float.h>
#include<math.h>
#include<string.h>
#include<ostream>

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
//and other functions to work properly and be able to throw errors
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



//dpois from dpois.c
double dpois_raw(double x, double lambda, int give_log)
// double dpois_raw(double x, double lambda, int log_p)
{
    /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
        lambda >= 0
    */
    if (lambda == 0) return( (x == 0) ? R_D__1 : R_D__0 );
    if (!R_FINITE(lambda)) return R_D__0; // including for the case where  x = lambda = +Inf
    if (x < 0) return( R_D__0 );
    if (x <= lambda * DBL_MIN) return(R_D_exp(-lambda) );
    if (lambda < x * DBL_MIN) {
	if (!R_FINITE(x)) // lambda < x = +Inf
	    return R_D__0;
	// else
	return(R_D_exp(-lambda + x*log(lambda) -lgammafn(x+1)));
    }
    return(R_D_fexp( M_2PI*x, -stirlerr(x)-bd0(x,lambda) ));
}



//from pgamma.c

/* dpois_wrap (x__1, lambda) := dpois(x__1 - 1, lambda);  where
 * dpois(k, L) := exp(-L) L^k / gamma(k+1)  {the usual Poisson probabilities}
 *
 * and  dpois*(.., give_log = TRUE) :=  log( dpois*(..) )
*/
static double
dpois_wrap (double x_plus_1, double lambda, int give_log)
{
#ifdef DEBUG_p
    REprintf (" dpois_wrap(x+1=%.14g, lambda=%.14g, log=%d)\n",
	      x_plus_1, lambda, give_log);
#endif
    if (!R_FINITE(lambda))
	return R_D__0;
    if (x_plus_1 > 1)
	return dpois_raw (x_plus_1 - 1, lambda, give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
	return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
	double d = dpois_raw (x_plus_1, lambda, give_log);
#ifdef DEBUG_p
	REprintf ("  -> d=dpois_raw(..)=%.14g\n", d);
#endif
	return give_log
	    ? d + log (x_plus_1 / lambda)
	    : d * (x_plus_1 / lambda);
    }
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


//pgamma_smallx from pgamma.c
/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
static double
pgamma_smallx (double x, double alph, int lower_tail, int log_p){
    double sum = 0, c = alph, n = 0, term;

#ifdef DEBUG_p
    REprintf (" pg_smallx(x=%.12g, alph=%.12g): ", x, alph);
#endif

    /*
     * Relative to 6.5.29 all terms have been multiplied by alph
     * and the first, thus being 1, is omitted.
     */

    do {
	n++;
	c *= -x / n;
	term = c / (alph + n);
	sum += term;
    } while (fabs (term) > DBL_EPSILON * fabs (sum));

#ifdef DEBUG_p
    REprintf ("%5.0f terms --> conv.sum=%g;", n, sum);
#endif
    if (lower_tail) {
	double f1 = log_p ? log1p (sum) : 1 + sum;
	double f2;
	if (alph > 1) {
	    f2 = dpois_raw (alph, x, log_p);
	    f2 = log_p ? f2 + x : f2 * exp (x);
	} else if (log_p)
	    f2 = alph * log (x) - lgamma1p (alph);
	else
	    f2 = pow (x, alph) / exp (lgamma1p (alph));
#ifdef DEBUG_p
    REprintf (" (f1,f2)= (%g,%g)\n", f1,f2);
#endif
	return log_p ? f1 + f2 : f1 * f2;
    } else {
	double lf2 = alph * log (x) - lgamma1p (alph);
#ifdef DEBUG_p
	REprintf (" 1:%.14g  2:%.14g\n", alph * log (x), lgamma1p (alph));
	REprintf (" sum=%.14g  log(1+sum)=%.14g	 lf2=%.14g\n",
		  sum, log1p (sum), lf2);
#endif
	if (log_p)
	    return R_Log1_Exp (log1p (sum) + lf2);
	else {
	    double f1m1 = sum;
	    double f2m1 = expm1 (lf2);
	    return -(f1m1 + f2m1 + f1m1 * f2m1);
	}
    }
} /* pgamma_smallx() */





//from pgamma.c
static double
pd_upper_series (double x, double y, int log_p){
    double term = x / y;
    double sum = term;

    do {
	y++;
	term *= x / y;
	sum += term;
    } while (term > sum * DBL_EPSILON);

    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
     *	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
     *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
     *	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
     */
    return log_p ? log (sum) : sum;
}


//from pgamma.c
/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
static double
pd_lower_cf (double y, double d){
    double f= 0.0 /* -Wall */, of, f0;
    double i, c2, c3, c4,  a1, b1,  a2, b2;

#define	NEEDED_SCALE				\
	  (b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

#ifdef DEBUG_p
    REprintf("pd_lower_cf(y=%.14g, d=%.14g)", y, d);
#endif
    if (y == 0) return 0;

    f0 = y/d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
    if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
#ifdef DEBUG_p
	REprintf(" very small 'y' -> returning (y/d)\n");
#endif
	return (f0);
    }

    if(f0 > 1.) f0 = 1.;
    c2 = y;
    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while NEEDED_SCALE

    i = 0; of = -1.; /* far away */
    while (i < max_it) {

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
	a1 = c4 * a2 + c3 * a1;
	b1 = c4 * b2 + c3 * b1;

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
	a2 = c4 * a1 + c3 * a2;
	b2 = c4 * b1 + c3 * b2;

	if NEEDED_SCALE

	if (b2 != 0) {
	    f = a2 / b2;
	    /* convergence check: relative; "absolute" for very small f : */
	    if (fabs (f - of) <= DBL_EPSILON * fmax2(f0, fabs(f))) {
#ifdef DEBUG_p
		REprintf(" %g iter.\n", i);
#endif
		return f;
	    }
	    of = f;
	}
    }

    MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
		    f);
    return f;/* should not happen ... */
} /* pd_lower_cf() */
#undef NEEDED_SCALE




//from pgamma.c
static double
pd_lower_series (double lambda, double y){
    double term = 1, sum = 0;

#ifdef DEBUG_p
    REprintf("pd_lower_series(lam=%.14g, y=%.14g) ...", lambda, y);
#endif
    while (y >= 1 && term > sum * DBL_EPSILON) {
	term *= y / lambda;
	sum += term;
	y--;
    }
    /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
     *	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
     *	   ~  y/lambda + o(y/lambda)
     */
#ifdef DEBUG_p
    REprintf(" done: term=%g, sum=%g, y= %g\n", term, sum, y);
#endif

    if (y != floor (y)) {
	/*
	 * The series does not converge as the terms start getting
	 * bigger (besides flipping sign) for y < -lambda.
	 */
	double f;
#ifdef DEBUG_p
	REprintf(" y not int: add another term ");
#endif
	/* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
	 *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
	f = pd_lower_cf (y, lambda + 1 - y);
#ifdef DEBUG_p
	REprintf("  (= %.14g) * term = %.14g to sum %g\n", f, term * f, sum);
#endif
	sum += term * f;
    }

    return sum;
} /* pd_lower_series() */




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




//ppois_asymp from pgamma.c
/*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */
static double
ppois_asymp (double x, double lambda, int lower_tail, int log_p){
    static const double coefs_a[8] = {
    -1e99, /* placeholder used for 1-indexing */
    2/3.,
    -4/135.,
    8/2835.,
    16/8505.,
    -8992/12629925.,
    -334144/492567075.,
    698752/1477701225.
    };

    static const double coefs_b[8] = {
    -1e99, /* placeholder */
    1/12.,
    1/288.,
    -139/51840.,
    -571/2488320.,
    163879/209018880.,
    5246819/75246796800.,
    -534703531/902961561600.
    };

    double elfb, elfb_term;
    double res12, res1_term, res1_ig, res2_term, res2_ig;
    double dfm, pt_, s2pt, f, np;
    int i;

    dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
    pt_ = - log1pmx (dfm / x);
    s2pt = sqrt (2 * x * pt_);
    if (dfm < 0) s2pt = -s2pt;

    res12 = 0;
    res1_ig = res1_term = sqrt (x);
    res2_ig = res2_term = s2pt;
    for (i = 1; i < 8; i++) {
    res12 += res1_ig * coefs_a[i];
    res12 += res2_ig * coefs_b[i];
    res1_term *= pt_ / i ;
    res2_term *= 2 * pt_ / (2 * i + 1);
    res1_ig = res1_ig / x + res1_term;
    res2_ig = res2_ig / x + res2_term;
    }

    elfb = x;
    elfb_term = 1;
    for (i = 1; i < 8; i++) {
    elfb += elfb_term * coefs_b[i];
    elfb_term /= x;
    }
    if (!lower_tail) elfb = -elfb;
#ifdef DEBUG_p
    REprintf ("res12 = %.14g   elfb=%.14g\n", elfb, res12);
#endif

    f = res12 / elfb;

    np = pnorm (s2pt, 0.0, 1.0, !lower_tail, log_p);

    if (log_p) {
    double n_d_over_p = dpnorm (s2pt, !lower_tail, np);
#ifdef DEBUG_p
    REprintf ("pp*_asymp(): f=%.14g	 np=e^%.14g  nd/np=%.14g  f*nd/np=%.14g\n",
    	  f, np, n_d_over_p, f * n_d_over_p);
#endif
    return np + log1p (f * n_d_over_p);
    } else {
    double nd = dnorm (s2pt, 0., 1., log_p);

#ifdef DEBUG_p
    REprintf ("pp*_asymp(): f=%.14g	 np=%.14g  nd=%.14g  f*nd=%.14g\n",
    	  f, np, nd, f * nd);
#endif
    return np + f * nd;
    }
} /* ppois_asymp() */






//////////////////////////////////////////////////////////////////////////////
//Now for pgamma functions
//////////////////////////////////////////////////////////////////////////////

double pgamma_raw (double x, double alph, int lower_tail, int log_p)
{
/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

    double res;

#ifdef DEBUG_p
    REprintf("pgamma_raw(x=%.14g, alph=%.14g, low=%d, log=%d)\n",
	     x, alph, lower_tail, log_p);
#endif
    R_P_bounds_01(x, 0., ML_POSINF);

    if (x < 1) {
	res = pgamma_smallx (x, alph, lower_tail, log_p);
    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
	/* incl. large alph compared to x */
	double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
	double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
	REprintf(" alph 'large': sum=pd_upper*()= %.12g, d=dpois_w(*)= %.12g\n",
		 sum, d);
#endif
	if (!lower_tail)
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
	else
	    res = log_p ? sum + d : sum * d;
    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
	/* incl. large x compared to alph */
	double sum;
	double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
	REprintf(" x 'large': d=dpois_w(*)= %.14g ", d);
#endif
	if (alph < 1) {
	    if (x * DBL_EPSILON > 1 - alph)
		sum = R_D__1;
	    else {
		double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
		/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
		sum = log_p ? log (f) : f;
	    }
	} else {
	    sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
	    sum = log_p ? log1p (sum) : 1 + sum;
	}
#ifdef DEBUG_p
	REprintf(", sum= %.14g\n", sum);
#endif
	if (!lower_tail)
	    res = log_p ? sum + d : sum * d;
	else
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
    } else { /* x >= 1 and x fairly near alph. */
#ifdef DEBUG_p
	REprintf(" using ppois_asymp()\n");
#endif
	res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
    }

    /*
     * We lose a fair amount of accuracy to underflow in the cases
     * where the final result is very close to DBL_MIN.	 In those
     * cases, simply redo via log space.
     */
    if (!log_p && res < DBL_MIN / DBL_EPSILON) {
	/* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
#ifdef DEBUG_p
	REprintf(" very small res=%.14g; -> recompute via log\n", res);
#endif
	return exp (pgamma_raw (x, alph, lower_tail, 1));
    } else
	return res;
}



double pgamma(double x, double alph, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(alph) || ISNAN(scale))
	return x + alph + scale;
#endif
    if(alph < 0. || scale <= 0.)
	ML_ERR_return_NAN;
    x /= scale;
#ifdef IEEE_754
    if (ISNAN(x)) /* eg. original x = scale = +Inf */
	return x;
#endif
    if(alph == 0.) /* limit case; useful e.g. in pnchisq() */
	return (x <= 0) ? R_DT_0: R_DT_1; /* <= assert  pgamma(0,0) ==> 0 */
    return pgamma_raw (x, alph, lower_tail, log_p);
}



//pchisq from pchisq.c
double pchisq(double x, double df, int lower_tail, int log_p)
{
    return pgamma(x, df/2., 2., lower_tail, log_p);
}

//////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////
//Elementary functions building up to pt
//////////////////////////////////////////////////////////////////////////////

int Rf_i1mach(int i)
{
    switch(i) {

    case  1: return 5;
    case  2: return 6;
    case  3: return 0;
    case  4: return 0;

    case  5: return CHAR_BIT * sizeof(int);
    case  6: return sizeof(int)/sizeof(char);

    case  7: return 2;
    case  8: return CHAR_BIT * sizeof(int) - 1;
    case  9: return INT_MAX;

    case 10: return FLT_RADIX;

    case 11: return FLT_MANT_DIG;
    case 12: return FLT_MIN_EXP;
    case 13: return FLT_MAX_EXP;

    case 14: return DBL_MANT_DIG;
    case 15: return DBL_MIN_EXP;
    case 16: return DBL_MAX_EXP;

    default: return 0;
    }
}



static double exparg(int l)
{
/* --------------------------------------------------------------------
 *     If l = 0 then  exparg(l) = The largest positive W for which
 *     exp(W) can be computed. With 0.99999 fuzz  ==> exparg(0) =   709.7756  nowadays

 *     if l = 1 (nonzero) then  exparg(l) = the largest negative W for
 *     which the computed value of exp(W) is nonzero.
 *     With 0.99999 fuzz			  ==> exparg(1) =  -709.0825  nowadays

 *     Note... only an approximate value for exparg(L) is needed.
 * -------------------------------------------------------------------- */

    static double const lnb = .69314718055995;
    int m = (l == 0) ? Rf_i1mach(16) : Rf_i1mach(15) - 1;

    return m * lnb * .99999;
} /* exparg */



//fpser from toms708.c
double fpser(double a, double b, double x, double eps, int log_p){
/* ----------------------------------------------------------------------- *

 *                 EVALUATION OF I (A,B)
 *                                X

 *          FOR B < MIN(EPS, EPS*A) AND X <= 0.5

 * ----------------------------------------------------------------------- */

    double ans, c, s, t, an, tol;

    /* SET  ans := x^a : */
    if (log_p) {
	ans = a * log(x);
    } else if (a > eps * 0.001) {
	t = a * log(x);
	if (t < exparg(1)) { /* exp(t) would underflow */
	    return 0.;
	}
	ans = exp(t);
    } else
	ans = 1.;

/*                NOTE THAT 1/B(A,B) = B */

    if (log_p)
	ans += log(b) - log(a);
    else
	ans *= b / a;

    tol = eps / a;
    an = a + 1.;
    t = x;
    s = t / an;
    do {
	an += 1.;
	t = x * t;
	c = t / an;
	s += c;
    } while (fabs(c) > tol);

    if (log_p)
	ans += log1p(a * s);
    else
	ans *= a * s + 1.;
    return ans;
} /* fpser */



static double psi(double x){
/* ---------------------------------------------------------------------

 *                 Evaluation of the Digamma function psi(x)

 *                           -----------

 *     Psi(xx) is assigned the value 0 when the digamma function cannot
 *     be computed.

 *     The main computation involves evaluation of rational Chebyshev
 *     approximations published in Math. Comp. 27, 123-127(1973) by
 *     Cody, Strecok and Thacher. */

/* --------------------------------------------------------------------- */
/*     Psi was written at Argonne National Laboratory for the FUNPACK */
/*     package of special function subroutines. Psi was modified by */
/*     A.H. Morris (NSWC). */
/* --------------------------------------------------------------------- */

    static double piov4 = .785398163397448; /* == pi / 4 */
/*     dx0 = zero of psi() to extended precision : */
    static double dx0 = 1.461632144968362341262659542325721325;

/* --------------------------------------------------------------------- */
/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) / (X - X0),  0.5 <= X <= 3. */
    static double p1[7] = { .0089538502298197,4.77762828042627,
	    142.441585084029,1186.45200713425,3633.51846806499,
	    4138.10161269013,1305.60269827897 };
    static double q1[6] = { 44.8452573429826,520.752771467162,
	    2210.0079924783,3641.27349079381,1908.310765963,
	    6.91091682714533e-6 };
/* --------------------------------------------------------------------- */


/* --------------------------------------------------------------------- */
/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) - LN(X) + 1 / (2*X),  X > 3. */

    static double p2[4] = { -2.12940445131011,-7.01677227766759,
	    -4.48616543918019,-.648157123766197 };
    static double q2[4] = { 32.2703493791143,89.2920700481861,
	    54.6117738103215,7.77788548522962 };
/* --------------------------------------------------------------------- */

    int i, m, n, nq;
    double d2;
    double w, z;
    double den, aug, sgn, xmx0, xmax11, upper, xsmall;

/* --------------------------------------------------------------------- */


/*     MACHINE DEPENDENT CONSTANTS ... */

/* --------------------------------------------------------------------- */
/*	  XMAX1	 = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
		   WITH ENTIRELY INT REPRESENTATION.  ALSO USED
		   AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
		   ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
		   PSI MAY BE REPRESENTED AS LOG(X).
 * Originally:  xmax11 = amin1(ipmpar(3), 1./spmpar(1))  */
    xmax11 = (double) INT_MAX;
    d2 = 0.5 / Rf_d1mach(3); /*= 0.5 / (0.5 * DBL_EPS) = 1/DBL_EPSILON = 2^52 */
    if(xmax11 > d2) xmax11 = d2;

/* --------------------------------------------------------------------- */
/*        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X) */
/*                 MAY BE REPRESENTED BY 1/X. */
    xsmall = 1e-9;
/* --------------------------------------------------------------------- */
    aug = 0.;
    if (x < 0.5) {
/* --------------------------------------------------------------------- */
/*     X < 0.5,  USE REFLECTION FORMULA */
/*     PSI(1-X) = PSI(X) + PI * COTAN(PI*X) */
/* --------------------------------------------------------------------- */
	if (fabs(x) <= xsmall) {

	    if (x == 0.) {
		goto L_err;
	    }
/* --------------------------------------------------------------------- */
/*     0 < |X| <= XSMALL.  USE 1/X AS A SUBSTITUTE */
/*     FOR  PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
	    aug = -1. / x;
	} else { /* |x| > xsmall */
/* --------------------------------------------------------------------- */
/*     REDUCTION OF ARGUMENT FOR COTAN */
/* --------------------------------------------------------------------- */
	    /* L100: */
	    w = -x;
	    sgn = piov4;
	    if (w <= 0.) {
		w = -w;
		sgn = -sgn;
	    }
/* --------------------------------------------------------------------- */
/*     MAKE AN ERROR EXIT IF |X| >= XMAX1 */
/* --------------------------------------------------------------------- */
	    if (w >= xmax11) {
		goto L_err;
	    }
	    nq = (int) w;
	    w -= (double) nq;
	    nq = (int) (w * 4.);
	    w = (w - (double) nq * 0.25) * 4.;
/* --------------------------------------------------------------------- */
/*     W IS NOW RELATED TO THE FRACTIONAL PART OF  4. * X. */
/*     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST */
/*     QUADRANT AND DETERMINE SIGN */
/* --------------------------------------------------------------------- */
	    n = nq / 2;
	    if (n + n != nq) {
		w = 1. - w;
	    }
	    z = piov4 * w;
	    m = n / 2;
	    if (m + m != n) {
		sgn = -sgn;
	    }
/* --------------------------------------------------------------------- */
/*     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
	    n = (nq + 1) / 2;
	    m = n / 2;
	    m += m;
	    if (m == n) {
/* --------------------------------------------------------------------- */
/*     CHECK FOR SINGULARITY */
/* --------------------------------------------------------------------- */
		if (z == 0.) {
		    goto L_err;
		}
/* --------------------------------------------------------------------- */
/*     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND */
/*     SIN/COS AS A SUBSTITUTE FOR TAN */
/* --------------------------------------------------------------------- */
		aug = sgn * (cos(z) / sin(z) * 4.);

	    } else { /* L140: */
		aug = sgn * (sin(z) / cos(z) * 4.);
	    }
	}

	x = 1. - x;

    }
    /* L200: */
    if (x <= 3.) {
/* --------------------------------------------------------------------- */
/*     0.5 <= X <= 3. */
/* --------------------------------------------------------------------- */
	den = x;
	upper = p1[0] * x;

	for (i = 1; i <= 5; ++i) {
	    den = (den + q1[i - 1]) * x;
	    upper = (upper + p1[i]) * x;
	}

	den = (upper + p1[6]) / (den + q1[5]);
	xmx0 = x - dx0;
	return den * xmx0 + aug;
    }

/* --------------------------------------------------------------------- */
/*     IF X >= XMAX1, PSI = LN(X) */
/* --------------------------------------------------------------------- */
    if (x < xmax11) {
/* --------------------------------------------------------------------- */
/*     3. < X < XMAX1 */
/* --------------------------------------------------------------------- */
	w = 1. / (x * x);
	den = w;
	upper = p2[0] * w;

	for (i = 1; i <= 3; ++i) {
	    den = (den + q2[i - 1]) * w;
	    upper = (upper + p2[i]) * w;
	}

	aug = upper / (den + q2[3]) - 0.5 / x + aug;
    }
    return aug + log(x);

/* --------------------------------------------------------------------- */
/*     ERROR RETURN */
/* --------------------------------------------------------------------- */
L_err:
    return 0.;
} /* psi */




static double apser(double a, double b, double x, double eps)
{
/* -----------------------------------------------------------------------
 *     apser() yields the incomplete beta ratio  I_{1-x}(b,a)  for
 *     a <= min(eps,eps*b), b*x <= 1, and x <= 0.5,  i.e., a is very small.
 *     Use only if above inequalities are satisfied.
 * ----------------------------------------------------------------------- */

    static double const g = .577215664901533;

    double tol, c, j, s, t, aj;
    double bx = b * x;

    t = x - bx;
    if (b * eps <= 0.02)
	c = log(x) + psi(b) + g + t;
    else // b > 2e13 : psi(b) ~= log(b)
	c = log(bx) + g + t;

    tol = eps * 5. * fabs(c);
    j = 1.;
    s = 0.;
    do {
	j += 1.;
	t *= x - bx / j;
	aj = t / j;
	s += aj;
    } while (fabs(aj) > tol);

    return -a * (c + s);
} /* apser */

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
//gam elementary functions from toms708.c
//////////////////////////////////////////////////////////////////////////////


//gam1 from toms708.c
static double gam1(double a)
{
/*     ------------------------------------------------------------------ */
/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5 */
/*     ------------------------------------------------------------------ */

    double d, t, w, bot, top;

    t = a;
    d = a - 0.5;
    // t := if(a > 1/2)  a-1  else  a
    if (d > 0.)
	t = d - 0.5;
    if (t < 0.) { /* L30: */
	static double
	    r[9] = { -.422784335098468,-.771330383816272,
		     -.244757765222226,.118378989872749,9.30357293360349e-4,
		     -.0118290993445146,.00223047661158249,2.66505979058923e-4,
		     -1.32674909766242e-4 },
	    s1 = .273076135303957,
	    s2 = .0559398236957378;

	top = (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]
		     ) * t + r[3]) * t + r[2]) * t + r[1]) * t + r[0];
	bot = (s2 * t + s1) * t + 1.;
	w = top / bot;
	R_ifDEBUG_printf("  gam1(a = %.15g): t < 0: w=%.15g\n", a, w);
	if (d > 0.)
	    return t * w / a;
	else
	    return a * (w + 0.5 + 0.5);

    } else if (t == 0) { // L10: a in {0, 1}
	return 0.;

    } else { /* t > 0;  L20: */
	static double
	    p[7] = { .577215664901533,-.409078193005776,
		     -.230975380857675,.0597275330452234,.0076696818164949,
		     -.00514889771323592,5.89597428611429e-4 },
	    q[5] = { 1.,.427569613095214,.158451672430138,
		     .0261132021441447,.00423244297896961 };

	top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]
		   ) * t + p[1]) * t + p[0];
	bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.;
	w = top / bot;
	R_ifDEBUG_printf("  gam1(a = %.15g): t > 0: (is a < 1.5 ?)  w=%.15g\n",
			 a, w);
	if (d > 0.) /* L21: */
	    return t / a * (w - 0.5 - 0.5);
	else
	    return a * w;
    }
} /* gam1 */






//gamln1 from toms708.c
static double gamln1(double a)
{
/* ----------------------------------------------------------------------- */
/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25 */
/* ----------------------------------------------------------------------- */

    double w;
    if (a < 0.6) {
	static double p0 = .577215664901533;
	static double p1 = .844203922187225;
	static double p2 = -.168860593646662;
	static double p3 = -.780427615533591;
	static double p4 = -.402055799310489;
	static double p5 = -.0673562214325671;
	static double p6 = -.00271935708322958;
	static double q1 = 2.88743195473681;
	static double q2 = 3.12755088914843;
	static double q3 = 1.56875193295039;
	static double q4 = .361951990101499;
	static double q5 = .0325038868253937;
	static double q6 = 6.67465618796164e-4;
	w = ((((((p6 * a + p5)* a + p4)* a + p3)* a + p2)* a + p1)* a + p0) /
	    ((((((q6 * a + q5)* a + q4)* a + q3)* a + q2)* a + q1)* a + 1.);
	return -(a) * w;
    }
    else { /* 0.6 <= a <= 1.25 */
	static double r0 = .422784335098467;
	static double r1 = .848044614534529;
	static double r2 = .565221050691933;
	static double r3 = .156513060486551;
	static double r4 = .017050248402265;
	static double r5 = 4.97958207639485e-4;
	static double s1 = 1.24313399877507;
	static double s2 = .548042109832463;
	static double s3 = .10155218743983;
	static double s4 = .00713309612391;
	static double s5 = 1.16165475989616e-4;
	double x = a - 0.5 - 0.5;
	w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
	    (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.);
	return x * w;
    }
} /* gamln1 */



//gamln from toms708.c
static double gamln(double a)
{
/* -----------------------------------------------------------------------
 *            Evaluation of  ln(gamma(a))  for positive a
 * ----------------------------------------------------------------------- */
/*     Written by Alfred H. Morris */
/*          Naval Surface Warfare Center */
/*          Dahlgren, Virginia */
/* ----------------------------------------------------------------------- */

    static double d = .418938533204673;/* d == 0.5*(LN(2*PI) - 1) */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    if (a <= 0.8)
	return gamln1(a) - log(a); /* ln(G(a+1)) - ln(a) == ln(G(a+1)/a) = ln(G(a)) */
    else if (a <= 2.25)
	return gamln1(a - 0.5 - 0.5);

    else if (a < 10.) {
	int i, n = (int)(a - 1.25);
	double t = a;
	double w = 1.;
	for (i = 1; i <= n; ++i) {
	    t += -1.;
	    w *= t;
	}
	return gamln1(t - 1.) + log(w);
    }
    else { /* a >= 10 */
	double t = 1. / (a * a);
	double w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
	return d + w + (a - 0.5) * (log(a) - 1.);
    }
} /* gamln */



//alnrel from toms708.c
static double alnrel(double a)
{
/* -----------------------------------------------------------------------
 *            Evaluation of the function ln(1 + a)
 * ----------------------------------------------------------------------- */

    if (fabs(a) > 0.375)
	return log(1. + a);
    // else : |a| <= 0.375
    static double
	p1 = -1.29418923021993,
	p2 = .405303492862024,
	p3 = -.0178874546012214,
	q1 = -1.62752256355323,
	q2 = .747811014037616,
	q3 = -.0845104217945565;
    double
	t = a / (a + 2.),
	t2 = t * t,
	w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.) /
	(((q3 * t2 + q2) * t2 + q1) * t2 + 1.);
    return t * 2. * w;

} /* alnrel */



//algdiv from toms708.c
static double algdiv(double a, double b)
{
/* ----------------------------------------------------------------------- */

/*     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8 */

/*                         -------- */

/*     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY */
/*     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X). */

/* ----------------------------------------------------------------------- */

    /* Initialized data */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    double c, d, h, t, u, v, w, x, s3, s5, x2, s7, s9, s11;

/* ------------------------ */
    if (a > b) {
	h = b / a;
	c = 1. / (h + 1.);
	x = h / (h + 1.);
	d = a + (b - 0.5);
    }
    else {
	h = a / b;
	c = h / (h + 1.);
	x = 1. / (h + 1.);
	d = b + (a - 0.5);
    }

/* Set s<n> = (1 - x^n)/(1 - x) : */

    x2 = x * x;
    s3 = x + x2 + 1.;
    s5 = x + x2 * s3 + 1.;
    s7 = x + x2 * s5 + 1.;
    s9 = x + x2 * s7 + 1.;
    s11 = x + x2 * s9 + 1.;

/* w := Del(b) - Del(a + b) */

    t = 1./ (b * b);
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
	    s3) * t + c0;
    w *= c / b;

/*                    COMBINE THE RESULTS */

    u = d * alnrel(a / b);
    v = a * (log(b) - 1.);
    if (u > v)
	return w - v - u;
    else
	return w - u - v;
} /* algdiv */



//bcorr from toms708.c
static double bcorr(double a0, double b0)
{
/* ----------------------------------------------------------------------- */

/*     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE */
/*     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A). */
/*     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8. */

/* ----------------------------------------------------------------------- */
    /* Initialized data */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    /* System generated locals */
    double ret_val, r1;

    /* Local variables */
    double a, b, c, h, t, w, x, s3, s5, x2, s7, s9, s11;
/* ------------------------ */
    a = min(a0, b0);
    b = max(a0, b0);

    h = a / b;
    c = h / (h + 1.);
    x = 1. / (h + 1.);
    x2 = x * x;

/*                SET SN = (1 - X^N)/(1 - X) */

    s3 = x + x2 + 1.;
    s5 = x + x2 * s3 + 1.;
    s7 = x + x2 * s5 + 1.;
    s9 = x + x2 * s7 + 1.;
    s11 = x + x2 * s9 + 1.;

/*                SET W = DEL(B) - DEL(A + B) */

/* Computing 2nd power */
    r1 = 1. / b;
    t = r1 * r1;
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
	    s3) * t + c0;
    w *= c / b;

/*                   COMPUTE  DEL(A) + W */

/* Computing 2nd power */
    r1 = 1. / a;
    t = r1 * r1;
    ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a +
	    w;
    return ret_val;
} /* bcorr */



//gsumln from toms708.c
static double gsumln(double a, double b)
{
/* ----------------------------------------------------------------------- */
/*          EVALUATION OF THE FUNCTION LN(GAMMA(A + B)) */
/*          FOR 1 <= A <= 2  AND  1 <= B <= 2 */
/* ----------------------------------------------------------------------- */

    double x = a + b - 2.;/* in [0, 2] */

    if (x <= 0.25)
	return gamln1(x + 1.);

    /* else */
    if (x <= 1.25)
	return gamln1(x) + alnrel(x);
    /* else x > 1.25 : */
    return gamln1(x - 1.) + log(x * (x + 1.));

} /* gsumln */




//betaln from toms708.c
static double betaln(double a0, double b0)
{
/* -----------------------------------------------------------------------
 *     Evaluation of the logarithm of the beta function  ln(beta(a0,b0))
 * ----------------------------------------------------------------------- */

    static double e = .918938533204673;/* e == 0.5*LN(2*PI) */

    double
	a = min(a0 ,b0),
	b = max(a0, b0);

    if (a < 8.) {
	if (a < 1.) {
/* ----------------------------------------------------------------------- */
//                    		A < 1
/* ----------------------------------------------------------------------- */
	    if (b < 8.)
		return gamln(a) + (gamln(b) - gamln(a+b));
	    else
		return gamln(a) + algdiv(a, b);
	}
	/* else */
/* ----------------------------------------------------------------------- */
//				1 <= A < 8
/* ----------------------------------------------------------------------- */
	double w;
    int n;
	if (a < 2.) {
	    if (b <= 2.) {
		return gamln(a) + gamln(b) - gsumln(a, b);
	    }
	    /* else */

	    if (b < 8.) {
		w = 0.;
		goto L40;
	    }
	    return gamln(a) + algdiv(a, b);
	}
	// else L30:    REDUCTION OF A WHEN B <= 1000

	if (b <= 1e3) {
	    n = (int)(a - 1.);
	    w = 1.;
	    for (int i = 1; i <= n; ++i) {
		a += -1.;
		double h = a / b;
		w *= h / (h + 1.);
	    }
	    w = log(w);

	    if (b >= 8.)
		return w + gamln(a) + algdiv(a, b);

	    // else
	L40:
	    // 	1 < A <= B < 8 :  reduction of B
	    n = (int)(b - 1.);
	    double z = 1.;
	    for (int i = 1; i <= n; ++i) {
		b += -1.;
		z *= b / (a + b);
	    }
	    return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)));
	}
	else { // L50:	reduction of A when  B > 1000
	    int n = (int)(a - 1.);
	    w = 1.;
	    for (int i = 1; i <= n; ++i) {
		a += -1.;
		w *= a / (a / b + 1.);
	    }
	    return log(w) - n * log(b) + (gamln(a) + algdiv(a, b));
	}

    } else {
/* ----------------------------------------------------------------------- */
	// L60:			A >= 8
/* ----------------------------------------------------------------------- */

	double
	    w = bcorr(a, b),
	    h = a / b,
	    u = -(a - 0.5) * log(h / (h + 1.)),
	    v = b * alnrel(h);
	if (u > v)
	    return log(b) * -0.5 + e + w - v - u;
	else
	    return log(b) * -0.5 + e + w - u - v;
    }

} /* betaln */



//bpser from toms708.c
static double bpser(double a, double b, double x, double eps, int log_p)
{
/* -----------------------------------------------------------------------
 * Power SERies expansion for evaluating I_x(a,b) when
 *	       b <= 1 or b*x <= 0.7.   eps is the tolerance used.
 * NB: if log_p is TRUE, also use it if   (b < 40  & lambda > 650)
 * ----------------------------------------------------------------------- */

    int i, m;
    double ans, c, t, u, z, a0, b0, apb;

    if (x == 0.) {
	return R_D__0;
    }
/* ----------------------------------------------------------------------- */
/*	      compute the factor  x^a/(a*Beta(a,b)) */
/* ----------------------------------------------------------------------- */
    a0 = min(a,b);
    if (a0 >= 1.) { /*		 ------	 1 <= a0 <= b0  ------ */
	z = a * log(x) - betaln(a, b);
	ans = log_p ? z - log(a) : exp(z) / a;
    }
    else {
	b0 = max(a,b);

	if (b0 < 8.) {

	    if (b0 <= 1.) { /*	 ------	 a0 < 1	 and  b0 <= 1  ------ */

		if(log_p) {
		    ans = a * log(x);
		} else {
		    ans = pow(x, a);
		    if (ans == 0.) /* once underflow, always underflow .. */
			return ans;
		}
		apb = a + b;
		if (apb > 1.) {
		    u = a + b - 1.;
		    z = (gam1(u) + 1.) / apb;
		} else {
		    z = gam1(apb) + 1.;
		}
		c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;

		if(log_p) /* FIXME ? -- improve quite a bit for c ~= 1 */
		    ans += log(c * (b / apb));
		else
		    ans *=  c * (b / apb);

	    } else { /* 	------	a0 < 1 < b0 < 8	 ------ */

		u = gamln1(a0);
		m = (int)(b0 - 1.);
		if (m >= 1) {
		    c = 1.;
		    for (i = 1; i <= m; ++i) {
			b0 += -1.;
			c *= b0 / (a0 + b0);
		    }
		    u += log(c);
		}

		z = a * log(x) - u;
		b0 += -1.; // => b0 in (0, 7)
		apb = a0 + b0;
		if (apb > 1.) {
		    u = a0 + b0 - 1.;
		    t = (gam1(u) + 1.) / apb;
		} else {
		    t = gam1(apb) + 1.;
		}

		if(log_p) /* FIXME? potential for improving log(t) */
		    ans = z + log(a0 / a) + log1p(gam1(b0)) - log(t);
		else
		    ans = exp(z) * (a0 / a) * (gam1(b0) + 1.) / t;
	    }

	} else { /* 		------  a0 < 1 < 8 <= b0  ------ */

	    u = gamln1(a0) + algdiv(a0, b0);
	    z = a * log(x) - u;

	    if(log_p)
		ans = z + log(a0 / a);
	    else
		ans = a0 / a * exp(z);
	}
    }
    R_ifDEBUG_printf(" bpser(a=%g, b=%g, x=%g, log=%d): prelim.ans = %.14g;\n",
		     a,b,x, log_p, ans);
    if (ans == R_D__0 || (!log_p && a <= eps * 0.1)) {
	return ans;
    }

/* ----------------------------------------------------------------------- */
/*		       COMPUTE THE SERIES */
/* ----------------------------------------------------------------------- */
    double tol = eps / a,
	n = 0.,
	sum = 0., w;
    c = 1.;
    do { // sum is alternating as long as n < b (<==> 1 - b/n < 0)
	n += 1.;
	c *= (0.5 - b / n + 0.5) * x;
	w = c / (a + n);
	sum += w;
    } while (n < 1e7 && fabs(w) > tol);
    if(fabs(w) > tol) { // the series did not converge (in time)
	// warn only when the result seems to matter:
	if(( log_p && !(a*sum > -1. && fabs(log1p(a * sum)) < eps*fabs(ans))) ||
	   (!log_p && fabs(a*sum + 1.) != 1.))
	    MATHLIB_WARNING5(
		" bpser(a=%g, b=%g, x=%g,...) did not converge (n=1e7, |w|/tol=%g > 1; A=%g)",
		a,b,x, fabs(w)/tol, ans);
    }
    R_ifDEBUG_printf("  -> n=%.0f iterations, |w|=%g %s %g=tol:=eps/a ==> a*sum=%g\n",
		     n, fabs(w), (fabs(w) > tol) ? ">!!>" : "<=",
		     tol, a*sum);
    if(log_p) {
	if (a*sum > -1.) ans += log1p(a * sum);
	else {
	    if(ans > ML_NEGINF)
		MATHLIB_WARNING3(
		    "pbeta(*, log.p=TRUE) -> bpser(a=%g, b=%g, x=%g,...) underflow to -Inf",
		    a,b,x);
	    ans = ML_NEGINF;
	}
    } else if (a*sum > -1.)
	ans *= (a * sum + 1.);
    else // underflow to
	ans = 0.;
    return ans;
} /* bpser */

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
//More elementary functions
//////////////////////////////////////////////////////////////////////////////


//esum in toms708.c
static double esum(int mu, double x, int give_log)
{
/* ----------------------------------------------------------------------- */
/*                    EVALUATION OF EXP(MU + X) */
/* ----------------------------------------------------------------------- */

    if(give_log)
	return x + (double) mu;

    // else :
    double w;
    if (x > 0.) { /* L10: */
	if (mu > 0)  return exp((double) mu) * exp(x);
	w = mu + x;
	if (w < 0.) return exp((double) mu) * exp(x);
    }
    else { /* x <= 0 */
	if (mu < 0)  return exp((double) mu) * exp(x);
	w = mu + x;
	if (w > 0.) return exp((double) mu) * exp(x);
    }
    return exp(w);

} /* esum */


//rlog1 in toms708.c
static double rlog1(double x)
{
/* -----------------------------------------------------------------------
 *             Evaluation of the function  x - ln(1 + x)
 * ----------------------------------------------------------------------- */

    static double a = .0566749439387324;
    static double b = .0456512608815524;
    static double p0 = .333333333333333;
    static double p1 = -.224696413112536;
    static double p2 = .00620886815375787;
    static double q1 = -1.27408923933623;
    static double q2 = .354508718369557;

    double h, r, t, w, w1;
    if (x < -0.39 || x > 0.57) { /* direct evaluation */
	w = x + 0.5 + 0.5;
	return x - log(w);
    }
    /* else */
    if (x < -0.18) { /* L10: */
	h = x + .3;
	h /= .7;
	w1 = a - h * .3;
    }
    else if (x > 0.18) { /* L20: */
	h = x * .75 - .25;
	w1 = b + h / 3.;
    }
    else { /*		Argument Reduction */
	h = x;
	w1 = 0.;
    }

/* L30:              	Series Expansion */

    r = h / (h + 2.);
    t = r * r;
    w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.);
    return t * 2. * (1. / (1. - r) - r * w) + w1;

} /* rlog1 */




//brcmp1 from toms708.c
//
// called only once from  bup(),  as   r = brcmp1(mu, a, b, x, y, FALSE) / a;
//                        -----
static double brcmp1(int mu, double a, double b, double x, double y, int give_log)
{
/* -----------------------------------------------------------------------
 *          Evaluation of    exp(mu) * x^a * y^b / beta(a,b)
 * ----------------------------------------------------------------------- */

    static double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
    /* R has  M_1_SQRT_2PI */

    /* Local variables */
    double c, t, u, v, z, a0, b0, apb;

    a0 = min(a,b);
    if (a0 < 8.) {
	double lnx, lny;
	if (x <= .375) {
	    lnx = log(x);
	    lny = alnrel(-x);
	} else if (y > .375) {
	    // L11:
	    lnx = log(x);
	    lny = log(y);
	} else {
	    lnx = alnrel(-y);
	    lny = log(y);
	}

	// L20:
	z = a * lnx + b * lny;
	if (a0 >= 1.) {
	    z -= betaln(a, b);
	    return esum(mu, z, give_log);
	}
	// else :
	/* ----------------------------------------------------------------------- */
	/*              PROCEDURE FOR A < 1 OR B < 1 */
	/* ----------------------------------------------------------------------- */
	// L30:
	b0 = max(a,b);
	if (b0 >= 8.) {
	/* L80:                  ALGORITHM FOR b0 >= 8 */
	    u = gamln1(a0) + algdiv(a0, b0);
	    R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1, b0 >= 8;  z=%.15g\n", z);
	    return give_log
		? log(a0) + esum(mu, z - u, true)
		:     a0  * esum(mu, z - u, false);

	} else if (b0 <= 1.) {
	    //                   a0 < 1, b0 <= 1
	    double ans = esum(mu, z, give_log);
	    if (ans == (give_log ? ML_NEGINF : 0.))
		return ans;

	    apb = a + b;
	    if (apb > 1.) {
		// L40:
		u = a + b - 1.;
		z = (gam1(u) + 1.) / apb;
	    } else {
		z = gam1(apb) + 1.;
	    }
	    // L50:
	    c = give_log
		? log1p(gam1(a)) + log1p(gam1(b)) - log(z)
		: (gam1(a) + 1.) * (gam1(b) + 1.) / z;
	    R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1, b0 <= 1;  c=%.15g\n", c);
	    return give_log
		? ans + log(a0) + c - log1p(a0 / b0)
		: ans * (a0 * c) / (a0 / b0 + 1.);
	}
	// else:               algorithm for	a0 < 1 < b0 < 8
	// L60:
	u = gamln1(a0);
	int n = (int)(b0 - 1.);
	if (n >= 1) {
	    c = 1.;
	    for (int i = 1; i <= n; ++i) {
		b0 += -1.;
		c *= b0 / (a0 + b0);
		/* L61: */
	    }
	    u += log(c); // TODO?: log(c) = log( prod(...) ) =  sum( log(...) )
	}
	// L70:
	z -= u;
	b0 += -1.;
	apb = a0 + b0;
	if (apb > 1.) {
	    // L71:
	    t = (gam1(apb - 1.) + 1.) / apb;
	} else {
	    t = gam1(apb) + 1.;
	}
	R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1 < b0 < 8;  t=%.15g\n", t);
	// L72:
	return give_log
	    ? log(a0)+ esum(mu, z, true) + log1p(gam1(b0)) - log(t) // TODO? log(t) = log1p(..)
	    :     a0 * esum(mu, z, false) * (gam1(b0) + 1.) / t;

    } else {

/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A >= 8 AND B >= 8 */
/* ----------------------------------------------------------------------- */
	// L100:
	double h, x0, y0, lambda;
	if (a > b) {
	    // L101:
	    h = b / a;
	    x0 = 1. / (h + 1.);// => lx0 := log(x0) = 0 - log1p(h)
	    y0 = h / (h + 1.);
	    lambda = (a + b) * y - b;
	} else {
	    h = a / b;
	    x0 = h / (h + 1.);  // => lx0 := log(x0) = - log1p(1/h)
	    y0 = 1. / (h + 1.);
	    lambda = a - (a + b) * x;
	}
	double lx0 = -log1p(b/a); // in both cases

	R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a,b >= 8;	x0=%.15g, lx0=log(x0)=%.15g\n",
			 x0, lx0);
	// L110:
	double e = -lambda / a;
	if (fabs(e) > 0.6) {
	    // L111:
	    u = e - log(x / x0);
	} else {
	    u = rlog1(e);
	}

	// L120:
	e = lambda / b;
	if (fabs(e) > 0.6) {
	    // L121:
	    v = e - log(y / y0);
	} else {
	    v = rlog1(e);
	}

	// L130:
	z = esum(mu, -(a * u + b * v), give_log);
	return give_log
	    ? log(const__)+ (log(b) + lx0)/2. + z      - bcorr(a, b)
	    :     const__ * sqrt(b * x0)      * z * exp(-bcorr(a, b));
    }

} /* brcmp1 */


//rexpm1 from toms708.c
double rexpm1(double x)
{
/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
/* ----------------------------------------------------------------------- */

    static double p1 = 9.14041914819518e-10;
    static double p2 = .0238082361044469;
    static double q1 = -.499999999085958;
    static double q2 = .107141568980644;
    static double q3 = -.0119041179760821;
    static double q4 = 5.95130811860248e-4;

    if (fabs(x) <= 0.15) {
	return x * (((p2 * x + p1) * x + 1.) /
		    ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.));
    }
    else { /* |x| > 0.15 : */
	double w = exp(x);
	if (x > 0.)
	    return w * (0.5 - 1. / w + 0.5);
	else
	    return w - 0.5 - 0.5;
    }

} /* rexpm1 */



//erf__ from toms708.c
static double erf__(double x)
{
/* -----------------------------------------------------------------------
 *             EVALUATION OF THE REAL ERROR FUNCTION
 * ----------------------------------------------------------------------- */

    /* Initialized data */

    static double c = .564189583547756;
    static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static double r[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* Local variables */
    double t, x2, ax, bot, top;

    ax = fabs(x);
    if (ax <= 0.5) {
	t = x * x;
	top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
	bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;

	return x * (top / bot);
    }

    // else:  |x| > 0.5

    if (ax <= 4.) { //  |x| in (0.5, 4]
	top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
		+ p[5]) * ax + p[6]) * ax + p[7];
	bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
		+ q[5]) * ax + q[6]) * ax + q[7];
	double R = 0.5 - exp(-x * x) * top / bot + 0.5;
	return (x < 0) ? -R : R;
    }

    // else:  |x| > 4

    if (ax >= 5.8) {
	return x > 0 ? 1 : -1;
    }

    // else:  4 < |x| < 5.8
    x2 = x * x;
    t = 1. / x2;
    top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
    t = (c - top / (x2 * bot)) / ax;
    double R = 0.5 - exp(-x2) * t + 0.5;
    return (x < 0) ? -R : R;
} /* erf */




//erfc1 from toms708.c
static double erfc1(int ind, double x)
{
/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
/* ----------------------------------------------------------------------- */

    /* Initialized data */

    static double c = .564189583547756;
    static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static double r[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    double ret_val;
    double e, t, w, bot, top;

    double ax = fabs(x);
    //				|X| <= 0.5 */
    if (ax <= 0.5) {
	double t = x * x,
	    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.,
	    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;
	ret_val = 0.5 - x * (top / bot) + 0.5;
	if (ind != 0) {
	    ret_val = exp(t) * ret_val;
	}
	return ret_val;
    }
    // else (L10:):		0.5 < |X| <= 4
    if (ax <= 4.) {
	top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
		+ p[5]) * ax + p[6]) * ax + p[7];
	bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
		+ q[5]) * ax + q[6]) * ax + q[7];
	ret_val = top / bot;

    } else { //			|X| > 4
	// L20:
	if (x <= -5.6) {
	    // L50:            	LIMIT VALUE FOR "LARGE" NEGATIVE X
	    ret_val = 2.;
	    if (ind != 0) {
		ret_val = exp(x * x) * 2.;
	    }
	    return ret_val;
	}
	if (ind == 0 && (x > 100. || x * x > -exparg(1))) {
	    // LIMIT VALUE FOR LARGE POSITIVE X   WHEN IND = 0
	    // L60:
	    return 0.;
	}

	// L30:
	t = 1. / (x * x);
	top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
	bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
	ret_val = (c - t * top / bot) / ax;
    }

    // L40:                 FINAL ASSEMBLY
    if (ind != 0) {
	if (x < 0.)
	    ret_val = exp(x * x) * 2. - ret_val;
    } else {
	// L41:  ind == 0 :
	w = x * x;
	t = w;
	e = w - t;
	ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
	if (x < 0.)
	    ret_val = 2. - ret_val;
    }
    return ret_val;

} /* erfc1 */



// called only from bgrat() , as   q_r = grat_r(b, z, log_r, eps)  :
static double grat_r(double a, double x, double log_r, double eps)
{
/* -----------------------------------------------------------------------
 *        Scaled complement of incomplete gamma ratio function
 *                   grat_r(a,x,r) :=  Q(a,x) / r
 * where
 *               Q(a,x) = pgamma(x,a, lower.tail=FALSE)
 *     and            r = e^(-x)* x^a / Gamma(a) ==  exp(log_r)
 *
 *     It is assumed that a <= 1.  eps is the tolerance to be used.
 * ----------------------------------------------------------------------- */

    if (a * x == 0.) { /* L130: */
	if (x <= a) {
	    /* L100: */ return exp(-log_r);
	} else {
	    /* L110:*/  return 0.;
	}
    }
    else if (a == 0.5) { // e.g. when called from pt()
	/* L120: */
	if (x < 0.25) {
	    double p = erf__(sqrt(x));
	    R_ifDEBUG_printf(" grat_r(a=%g, x=%g ..)): a=1/2 --> p=erf__(.)= %g\n",
			     a, x, p);
	    return (0.5 - p + 0.5)*exp(-log_r);

        } else { // 2013-02-27: improvement for "large" x: direct computation of q/r:
	    double sx = sqrt(x),
		q_r = erfc1(1, sx)/sx * M_SQRT_PI;
	    R_ifDEBUG_printf(" grat_r(a=%g, x=%g ..)): a=1/2 --> q_r=erfc1(..)/r= %g\n",
			     a,x, q_r);
	    return q_r;
	}

    } else if (x < 1.1) { /* L10:  Taylor series for  P(a,x)/x^a */

	double an = 3.,
	    c = x,
	    sum = x / (a + 3.),
	    tol = eps * 0.1 / (a + 1.), t;
	do {
	    an += 1.;
	    c *= -(x / an);
	    t = c / (a + an);
	    sum += t;
	} while (fabs(t) > tol);

	R_ifDEBUG_printf(" grat_r(a=%g, x=%g, log_r=%g): sum=%g; Taylor w/ %.0f terms",
			 a,x,log_r, sum, an-3.);
	double j = a * x * ((sum/6. - 0.5/(a + 2.)) * x + 1./(a + 1.)),
	    z = a * log(x),
	    h = gam1(a),
	    g = h + 1.;

	if ((x >= 0.25 && (a < x / 2.59)) || (z > -0.13394)) {
	    // L40:
	    double l = rexpm1(z),
		q = ((l + 0.5 + 0.5) * j - l) * g - h;
	    if (q <= 0.) {
		R_ifDEBUG_printf(" => q_r= 0.\n");
		/* L110:*/ return 0.;
	    } else {
		R_ifDEBUG_printf(" => q_r=%.15g\n", q * exp(-log_r));
		return q * exp(-log_r);
	    }

	} else {
	    double p = exp(z) * g * (0.5 - j + 0.5);
	    R_ifDEBUG_printf(" => q_r=%.15g\n", (0.5 - p + 0.5) * exp(-log_r));
	    return /* q/r = */ (0.5 - p + 0.5) * exp(-log_r);
	}

    } else {
	/* L50: ----  (x >= 1.1)  ---- Continued Fraction Expansion */

	double a2n_1 = 1.,
	    a2n = 1.,
	    b2n_1 = x,
	    b2n = x + (1. - a),
	    c = 1., am0, an0;

	do {
	    a2n_1 = x * a2n + c * a2n_1;
	    b2n_1 = x * b2n + c * b2n_1;
	    am0 = a2n_1 / b2n_1;
	    c += 1.;
	    double c_a = c - a;
	    a2n = a2n_1 + c_a * a2n;
	    b2n = b2n_1 + c_a * b2n;
	    an0 = a2n / b2n;
	} while (fabs(an0 - am0) >= eps * an0);

	R_ifDEBUG_printf(" grat_r(a=%g, x=%g, log_r=%g): Cont.frac. %.0f terms => q_r=%.15g\n",
			 a,x, log_r, c-1., an0);
	return /* q/r = (r * an0)/r = */ an0;
    }
} /* grat_r */



//bup from toms708.c
static double bup(double a, double b, double x, double y, int n, double eps,
		  int give_log)
{
/* ----------------------------------------------------------------------- */
/*     EVALUATION OF I_x(A,B) - I_x(A+N,B) WHERE N IS A POSITIVE INT. */
/*     EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */

    double ret_val;
    int i, k, mu;
    double d, l;

// Obtain the scaling factor exp(-mu) and exp(mu)*(x^a * y^b / beta(a,b))/a

    double apb = a + b,
	ap1 = a + 1.;
    if (n > 1 && a >= 1. && apb >= ap1 * 1.1) {
	mu = (int)fabs(exparg(1));
	k = (int) exparg(0);
	if (mu > k)
	    mu = k;
 	d = exp(-(double) mu);
    }
    else {
	mu = 0;
	d = 1.;
    }

    /* L10: */
    ret_val = give_log
	? brcmp1(mu, a, b, x, y, true) - log(a)
	: brcmp1(mu, a, b, x, y, false)  / a;
    if (n == 1 ||
	(give_log && ret_val == ML_NEGINF) || (!give_log && ret_val == 0.))
	return ret_val;

    int nm1 = n - 1;
    double w = d;

/*          LET K BE THE INDEX OF THE MAXIMUM TERM */

    k = 0;
    if (b > 1.) {
	if (y > 1e-4) {
	    double r = (b - 1.) * x / y - a;
	    if (r >= 1.)
		k = (r < nm1) ? (int) r : nm1;
	} else
	    k = nm1;

//          ADD THE INCREASING TERMS OF THE SERIES - if k > 0
/* L30: */
	for (i = 0; i < k; ++i) {
	    l = (double) i;
	    d *= (apb + l) / (ap1 + l) * x;
	    w += d;
	}
    }

// L40:     ADD THE REMAINING TERMS OF THE SERIES

    for (i = k; i < nm1; ++i) {
	l = (double) i;
	d *= (apb + l) / (ap1 + l) * x;
	w += d;
	if (d <= eps * w) /* relativ convergence (eps) */
	    break;
    }

    // L50: TERMINATE THE PROCEDURE
    if(give_log) {
	ret_val += log(w);
    } else
	ret_val *= w;

    return ret_val;
} /* bup */


//logspace_add from pgamma.c
/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_add (double logx, double logy)
{
    return fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
}


//bgrat from toms708.c
static void bgrat(double a, double b, double x, double y, double *w,
		  double eps, int *ierr, Rboolean log_w)
{
/* -----------------------------------------------------------------------
*     Asymptotic Expansion for I_x(a,b)  when a is larger than b.
*     Compute   w := w + I_x(a,b)
*     It is assumed a >= 15 and b <= 1.
*     eps is the tolerance used.
*     ierr is a variable that reports the status of the results.
*
* if(log_w),  *w  itself must be in log-space;
*     compute   w := w + I_x(a,b)  but return *w = log(w):
*          *w := log(exp(*w) + I_x(a,b)) = logspace_add(*w, log( I_x(a,b) ))
* ----------------------------------------------------------------------- */

#define n_terms_bgrat 30
    double c[n_terms_bgrat], d[n_terms_bgrat];
    double bm1 = b - 0.5 - 0.5,
	nu = a + bm1 * 0.5, /* nu = a + (b-1)/2 =: T, in (9.1) of
			     * Didonato & Morris(1992), p.362 */
	lnx = (y > 0.375) ? log(x) : alnrel(-y),
	z = -nu * lnx; // z =: u in (9.1) of D.&M.(1992)

    if (b * z == 0.) { // should not happen, but does, e.g.,
	// for  pbeta(1e-320, 1e-5, 0.5)  i.e., _subnormal_ x,
	// Warning ... bgrat(a=20.5, b=1e-05, x=1, y=9.99989e-321): ..
	MATHLIB_WARNING5(
	    "bgrat(a=%g, b=%g, x=%g, y=%g): z=%g, b*z == 0 underflow, hence inaccurate pbeta()",
	    a,b,x,y, z);
	/* L_Error:    THE EXPANSION CANNOT BE COMPUTED */
	 *ierr = 1; return;
    }

/*                 COMPUTATION OF THE EXPANSION */
    double
	/* r1 = b * (gam1(b) + 1.) * exp(b * log(z)),// = b/gamma(b+1) z^b = z^b / gamma(b)
	 * set r := exp(-z) * z^b / gamma(b) ;
	 *          gam1(b) = 1/gamma(b+1) - 1 , b in [-1/2, 3/2] */
	// exp(a*lnx) underflows for large (a * lnx); e.g. large a ==> using log_r := log(r):
	// r = r1 * exp(a * lnx) * exp(bm1 * 0.5 * lnx);
	// log(r)=log(b) + log1p(gam1(b)) + b * log(z) + (a * lnx) + (bm1 * 0.5 * lnx),
	log_r = log(b) + log1p(gam1(b)) + b * log(z) + nu * lnx,
	// FIXME work with  log_u = log(u)  also when log_p=FALSE  (??)
	// u is 'factored out' from the expansion {and multiplied back, at the end}:
	log_u = log_r - (algdiv(b, a) + b * log(nu)),// algdiv(b,a) = log(gamma(a)/gamma(a+b))
	/* u = (log_p) ? log_r - u : exp(log_r-u); // =: M  in (9.2) of {reference above} */
	/* u = algdiv(b, a) + b * log(nu);// algdiv(b,a) = log(gamma(a)/gamma(a+b)) */
	// u = (log_p) ? log_u : exp(log_u); // =: M  in (9.2) of {reference above}
	u = exp(log_u);

    if (log_u == ML_NEGINF) {
	R_ifDEBUG_printf(" bgrat(*): underflow log_u = -Inf  = log_r -u', log_r = %g ",
			 log_r);
	/* L_Error:    THE EXPANSION CANNOT BE COMPUTED */ *ierr = 2; return;
    }

    Rboolean u_0 = (u == 0.); // underflow --> do work with log(u) == log_u !
    double l = // := *w/u .. but with care: such that it also works when u underflows to 0:
	log_w
	? ((*w == ML_NEGINF) ? 0. : exp(  *w    - log_u))
	: ((*w == 0.)        ? 0. : exp(log(*w) - log_u));

    R_ifDEBUG_printf(" bgrat(a=%g, b=%g, x=%g, *)\n -> u=%g, l='w/u'=%g, ",
		     a,b,x, u, l);
    double
	q_r = grat_r(b, z, log_r, eps), // = q/r of former grat1(b,z, r, &p, &q)
	v = 0.25 / (nu * nu),
	t2 = lnx * 0.25 * lnx,
	j = q_r,
	sum = j,
	t = 1., cn = 1., n2 = 0.;
    for (int n = 1; n <= n_terms_bgrat; ++n) {
	double bp2n = b + n2;
	j = (bp2n * (bp2n + 1.) * j + (z + bp2n + 1.) * t) * v;
	n2 += 2.;
	t *= t2;
	cn /= n2 * (n2 + 1.);
	int nm1 = n - 1;
	c[nm1] = cn;
	double s = 0.;
	if (n > 1) {
	    double coef = b - n;
	    for (int i = 1; i <= nm1; ++i) {
		s += coef * c[i - 1] * d[nm1 - i];
		coef += b;
	    }
	}
	d[nm1] = bm1 * cn + s / n;
	double dj = d[nm1] * j;
	sum += dj;
	if (sum <= 0.) {
	    R_ifDEBUG_printf(" bgrat(*): sum_n(..) <= 0; should not happen (n=%d)\n", n);
	    /* L_Error:    THE EXPANSION CANNOT BE COMPUTED */ *ierr = 3; return;
	}
	if (fabs(dj) <= eps * (sum + l)) {
	    *ierr = 0;
	    break;
	} else if(n == n_terms_bgrat) { // never? ; please notify R-core if seen:
	    *ierr = 4;
	    MATHLIB_WARNING5(
	"bgrat(a=%g, b=%g, x=%g) *no* convergence: NOTIFY R-core!\n dj=%g, rel.err=%g\n",
		a,b,x, dj, fabs(dj) /(sum + l));
	}
    } // for(n .. n_terms..)

/*                    ADD THE RESULTS TO W */

    if(log_w) // *w is in log space already:
	*w = logspace_add(*w, log_u + log(sum));
    else
	*w += (u_0 ? exp(log_u + log(sum)) : u * sum);
    return;
} /* bgrat */




//basym from toms708.c
static double basym(double a, double b, double lambda, double eps, int log_p)
{
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR I_x(A,B) FOR LARGE A AND B. */
/*     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED. */
/*     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT */
/*     A AND B ARE GREATER THAN OR EQUAL TO 15. */
/* ----------------------------------------------------------------------- */


/* ------------------------ */
/*     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
/*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
#define num_IT 20
/*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */

    static double const e0 = 1.12837916709551;/* e0 == 2/sqrt(pi) */
    static double const e1 = .353553390593274;/* e1 == 2^(-3/2)   */
    static double const ln_e0 = 0.120782237635245; /* == ln(e0) */

    double a0[num_IT + 1], b0[num_IT + 1], c[num_IT + 1], d[num_IT + 1];

    double f = a * rlog1(-lambda/a) + b * rlog1(lambda/b), t;
    if(log_p)
	t = -f;
    else {
	t = exp(-f);
	if (t == 0.) {
	    return 0; /* once underflow, always underflow .. */
	}
    }
    double z0 = sqrt(f),
	z = z0 / e1 * 0.5,
	z2 = f + f,
	h, r0, r1, w0;

    if (a < b) {
	h = a / b;
	r0 = 1. / (h + 1.);
	r1 = (b - a) / b;
	w0 = 1. / sqrt(a * (h + 1.));
    } else {
	h = b / a;
	r0 = 1. / (h + 1.);
	r1 = (b - a) / a;
	w0 = 1. / sqrt(b * (h + 1.));
    }

    a0[0] = r1 * .66666666666666663;
    c[0] = a0[0] * -0.5;
    d[0] = -c[0];
    double j0 = 0.5 / e0 * erfc1(1, z0),
	j1 = e1,
	sum = j0 + d[0] * w0 * j1;

    double s = 1.,
	h2 = h * h,
	hn = 1.,
	w = w0,
	znm1 = z,
	zn = z2;
    for (int n = 2; n <= num_IT; n += 2) {
	hn *= h2;
	a0[n - 1] = r0 * 2. * (h * hn + 1.) / (n + 2.);
	int np1 = n + 1;
	s += hn;
	a0[np1 - 1] = r1 * 2. * s / (n + 3.);

	for (int i = n; i <= np1; ++i) {
	    double r = (i + 1.) * -0.5;
	    b0[0] = r * a0[0];
	    for (int m = 2; m <= i; ++m) {
		double bsum = 0.;
		for (int j = 1; j <= m-1; ++j) {
		    int mmj = m - j;
		    bsum += (j * r - mmj) * a0[j - 1] * b0[mmj - 1];
		}
		b0[m - 1] = r * a0[m - 1] + bsum / m;
	    }
	    c[i - 1] = b0[i - 1] / (i + 1.);

	    double dsum = 0.;
	    for (int j = 1; j <= i-1; ++j) {
		dsum += d[i - j - 1] * c[j - 1];
	    }
	    d[i - 1] = -(dsum + c[i - 1]);
	}

	j0 = e1 * znm1 + (n - 1.) * j0;
	j1 = e1 * zn + n * j1;
	znm1 = z2 * znm1;
	zn = z2 * zn;
	w *= w0;
	double t0 = d[n - 1] * w * j0;
	w *= w0;
	double t1 = d[np1 - 1] * w * j1;
	sum += t0 + t1;
	if (fabs(t0) + fabs(t1) <= eps * sum) {
	    break;
	}
    }

    if(log_p)
	return ln_e0 + t - bcorr(a, b) + log(sum);
    else {
	double u = exp(-bcorr(a, b));
	return e0 * t * u * sum;
    }

} /* basym_ */


//brcomp from toms708.c
static double brcomp(double a, double b, double x, double y, int log_p)
{
/* -----------------------------------------------------------------------
 *		 Evaluation of x^a * y^b / Beta(a,b)
 * ----------------------------------------------------------------------- */

    static double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
    /* R has  M_1_SQRT_2PI , and M_LN_SQRT_2PI = ln(sqrt(2*pi)) = 0.918938.. */
    int i, n;
    double c, e, u, v, z, a0, b0, apb;

    if (x == 0. || y == 0.) {
	return R_D__0;
    }
    a0 = min(a, b);
    if (a0 < 8.) {
	double lnx, lny;
	if (x <= .375) {
	    lnx = log(x);
	    lny = alnrel(-x);
	}
	else {
	    if (y > .375) {
		lnx = log(x);
		lny = log(y);
	    } else {
		lnx = alnrel(-y);
		lny = log(y);
	    }
	}

	z = a * lnx + b * lny;
	if (a0 >= 1.) {
	    z -= betaln(a, b);
	    return R_D_exp(z);
	}

/* ----------------------------------------------------------------------- */
/*		PROCEDURE FOR a < 1 OR b < 1 */
/* ----------------------------------------------------------------------- */

	b0 = max(a, b);
	if (b0 >= 8.) { /* L80: */
	    u = gamln1(a0) + algdiv(a0, b0);

	    return (log_p ? log(a0) + (z - u)  : a0 * exp(z - u));
	}
	/* else : */

	if (b0 <= 1.) { /*		algorithm for max(a,b) = b0 <= 1 */

	    double e_z = R_D_exp(z);

	    if (!log_p && e_z == 0.) /* exp() underflow */
		return 0.;

	    apb = a + b;
	    if (apb > 1.) {
		u = a + b - 1.;
		z = (gam1(u) + 1.) / apb;
	    } else {
		z = gam1(apb) + 1.;
	    }

	    c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;
	    /* FIXME? log(a0*c)= log(a0)+ log(c) and that is improvable */
	    return (log_p
		    ? e_z + log(a0 * c) - log1p(a0/b0)
		    : e_z * (a0 * c) / (a0 / b0 + 1.));
	}

	/* else : 		  ALGORITHM FOR 1 < b0 < 8 */

	u = gamln1(a0);
	n = (int)(b0 - 1.);
	if (n >= 1) {
	    c = 1.;
	    for (i = 1; i <= n; ++i) {
		b0 += -1.;
		c *= b0 / (a0 + b0);
	    }
	    u = log(c) + u;
	}
	z -= u;
	b0 += -1.;
	apb = a0 + b0;
	double t;
	if (apb > 1.) {
	    u = a0 + b0 - 1.;
	    t = (gam1(u) + 1.) / apb;
	} else {
	    t = gam1(apb) + 1.;
	}

	return (log_p
		? log(a0) + z + log1p(gam1(b0))  - log(t)
		: a0 * exp(z) * (gam1(b0) + 1.) / t);

    } else {
/* ----------------------------------------------------------------------- */
/*		PROCEDURE FOR A >= 8 AND B >= 8 */
/* ----------------------------------------------------------------------- */
	double h, x0, y0, lambda;
	if (a <= b) {
	    h = a / b;
	    x0 = h / (h + 1.);
	    y0 = 1. / (h + 1.);
	    lambda = a - (a + b) * x;
	} else {
	    h = b / a;
	    x0 = 1. / (h + 1.);
	    y0 = h / (h + 1.);
	    lambda = (a + b) * y - b;
	}

	e = -lambda / a;
	if (fabs(e) > .6)
	    u = e - log(x / x0);
	else
	    u = rlog1(e);

	e = lambda / b;
	if (fabs(e) <= .6)
	    v = rlog1(e);
	else
	    v = e - log(y / y0);

	z = log_p ? -(a * u + b * v) : exp(-(a * u + b * v));

	return(log_p
	       ? -M_LN_SQRT_2PI + .5*log(b * x0) + z - bcorr(a,b)
	       : const__ * sqrt(b * x0) * z * exp(-bcorr(a, b)));
    }
} /* brcomp */



//bfrac from toms708.c
static double bfrac(double a, double b, double x, double y, double lambda,
		    double eps, int log_p)
{
/* -----------------------------------------------------------------------
       Continued fraction expansion for I_x(a,b) when a, b > 1.
       It is assumed that  lambda = (a + b)*y - b.
   -----------------------------------------------------------------------*/

    double c, e, n, p, r, s, t, w, c0, c1, r0, an, bn, yp1, anp1, bnp1,
	beta, alpha, brc;

    if(!R_FINITE(lambda)) return ML_NAN;// TODO: can return 0 or 1 (?)
    R_ifDEBUG_printf(" bfrac(a=%g, b=%g, x=%g, y=%g, lambda=%g, eps=%g, log_p=%d):",
		     a,b,x,y, lambda, eps, log_p);
    brc = brcomp(a, b, x, y, log_p);
    if(ISNAN(brc)) { // e.g. from   L <- 1e308; pnbinom(L, L, mu = 5)
	R_ifDEBUG_printf(" --> brcomp(a,b,x,y) = NaN\n");
	ML_ERR_return_NAN; // TODO: could we know better?
    }
    if (!log_p && brc == 0.) {
	R_ifDEBUG_printf(" --> brcomp(a,b,x,y) underflowed to 0.\n");
	return 0.;
    }
#ifdef DEBUG_bratio
    else
	REprintf("\n");
#endif

    c = lambda + 1.;
    c0 = b / a;
    c1 = 1. / a + 1.;
    yp1 = y + 1.;

    n = 0.;
    p = 1.;
    s = a + 1.;
    an = 0.;
    bn = 1.;
    anp1 = 1.;
    bnp1 = c / c1;
    r = c1 / c;

/*        CONTINUED FRACTION CALCULATION */

    do {
	n += 1.;
	t = n / a;
	w = n * (b - n) * x;
	e = a / s;
	alpha = p * (p + c0) * e * e * (w * x);
	e = (t + 1.) / (c1 + t + t);
	beta = n + w / s + e * (c + n * yp1);
	p = t + 1.;
	s += 2.;

	/* update an, bn, anp1, and bnp1 */

	t = alpha * an + beta * anp1;	an = anp1;	anp1 = t;
	t = alpha * bn + beta * bnp1;	bn = bnp1;	bnp1 = t;

	r0 = r;
	r = anp1 / bnp1;
#ifdef _not_normally_DEBUG_bfrac
	R_ifDEBUG_printf(" n=%5.0f, a_{n,n+1}= (%12g,%12g),  b_{n,n+1} = (%12g,%12g) => r0,r = (%14g,%14g)\n",
			 n, an,anp1, bn,bnp1, r0, r);
#endif
	if (fabs(r - r0) <= eps * r)
	    break;

	/* rescale an, bn, anp1, and bnp1 */

	an /= bnp1;
	bn /= bnp1;
	anp1 = r;
	bnp1 = 1.;
    } while (n < 10000);// arbitrary; had '1' --> infinite loop for  lambda = Inf
    R_ifDEBUG_printf("  in bfrac(): n=%.0f terms cont.frac.; brc=%g, r=%g\n",
		     n, brc, r);
    if(n >= 10000 && fabs(r - r0) > eps * r)
	MATHLIB_WARNING5(
	    " bfrac(a=%g, b=%g, x=%g, y=%g, lambda=%g) did *not* converge (in 10000 steps)\n",
	    a,b,x,y, lambda);
    return (log_p ? brc + log(r) : brc * r);
} /* bfrac */


void 
bratio(double a, double b, double x, double y, double *w, double *w1,
       int *ierr, int log_p)
{

    Rboolean do_swap;
    int n, ierr1 = 0;
    double z, a0, b0, x0, y0, lambda;

/*  eps is a machine dependent constant: the smallest
 *      floating point number for which   1. + eps > 1.
 * NOTE: for almost all purposes it is replaced by 1e-15 (~= 4.5 times larger) below */
    double eps = 2. * Rf_d1mach(3); /* == DBL_EPSILON (in R, Rmath) */

/* ----------------------------------------------------------------------- */
    *w  = R_D__0;
    *w1 = R_D__0;

#ifdef IEEE_754
    // safeguard, preventing infinite loops further down
    if (ISNAN(x) || ISNAN(y) ||
	ISNAN(a) || ISNAN(b)) { *ierr = 9; return; }
#endif
    if (a < 0. || b < 0.)   { *ierr = 1; return; }
    if (a == 0. && b == 0.) { *ierr = 2; return; }
    if (x < 0. || x > 1.)   { *ierr = 3; return; }
    if (y < 0. || y > 1.)   { *ierr = 4; return; }

    /* check that  'y == 1 - x' : */
    z = x + y - 0.5 - 0.5;

    if (fabs(z) > eps * 3.) { *ierr = 5; return; }

    R_ifDEBUG_printf("bratio(a=%g, b=%g, x=%9g, y=%9g, .., log_p=%d): ",
		     a,b,x,y, log_p);
    *ierr = 0;

    Rboolean a_lt_b = (a < b); //line was moved

    if (x == 0.) goto L200;
    if (y == 0.) goto L210;

    if (a == 0.) goto L211;
    if (b == 0.) goto L201;
    eps = max(eps, 1e-15);
    
    if (/* max(a,b) */ (a_lt_b ? b : a) < eps * .001) { /* procedure for a and b < 0.001 * eps */
	// L230:  -- result *independent* of x (!)
	// *w  = a/(a+b)  and  w1 = b/(a+b) :
	if(log_p) {
	    if(a_lt_b) {
		*w  = log1p(-a/(a+b)); // notably if a << b
		*w1 = log  ( a/(a+b));
	    } else { // b <= a
		*w  = log  ( b/(a+b));
		*w1 = log1p(-b/(a+b));
	    }
	} else {
	    *w	= b / (a + b);
	    *w1 = a / (a + b);
	}

	R_ifDEBUG_printf("a & b very small -> simple ratios (%g,%g)\n", *w,*w1);
	return;
    }

#define SET_0_noswap \
    a0 = a;  x0 = x; \
    b0 = b;  y0 = y;

#define SET_0_swap   \
    a0 = b;  x0 = y; \
    b0 = a;  y0 = x;

    if (min(a,b) <= 1.) { /*------------------------ a <= 1  or  b <= 1 ---- */

	do_swap = (x > 0.5);
	if (do_swap) {
	    SET_0_swap;
	} else {
	    SET_0_noswap;
	}
	/* now have  x0 <= 1/2 <= y0  (still  x0+y0 == 1) */

	R_ifDEBUG_printf(" min(a,b) <= 1, do_swap=%d;", do_swap);

	if (b0 < min(eps, eps * a0)) { /* L80: */
	    *w = fpser(a0, b0, x0, eps, log_p);
	    *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
	    R_ifDEBUG_printf("  b0 small -> w := fpser(*) = %.15g\n", *w);
	    goto L_end;
	}

	if (a0 < min(eps, eps * b0) && b0 * x0 <= 1.) { /* L90: */
	    *w1 = apser(a0, b0, x0, eps);
	    R_ifDEBUG_printf("  a0 small -> w1 := apser(*) = %.15g\n", *w1);
	    goto L_end_from_w1;
	}

	Rboolean did_bup = false;
	if (max(a0,b0) > 1.) { /* L20:  min(a,b) <= 1 < max(a,b)  */
	    R_ifDEBUG_printf("\n L20:  min(a,b) <= 1 < max(a,b); ");
	    if (b0 <= 1.) goto L_w_bpser;

	    if (x0 >= 0.29) /* was 0.3, PR#13786 */	goto L_w1_bpser;

	    if (x0 < 0.1 && pow(x0*b0, a0) <= 0.7)	goto L_w_bpser;

	    if (b0 > 15.) {
		*w1 = 0.;
		goto L131;
	    }
	} else { /*  a, b <= 1 */
	    R_ifDEBUG_printf("\n      both a,b <= 1; ");
	    if (a0 >= min(0.2, b0))	goto L_w_bpser;

	    if (pow(x0, a0) <= 0.9) 	goto L_w_bpser;

	    if (x0 >= 0.3)		goto L_w1_bpser;
	}
	n = 20; /* goto L130; */
	*w1 = bup(b0, a0, y0, x0, n, eps, false); did_bup = true;
	R_ifDEBUG_printf("  ... n=20 and *w1 := bup(*) = %.15g; ", *w1);
	b0 += n;
    L131:
	R_ifDEBUG_printf(" L131: bgrat(*, w1=%.15g) ", *w1);
	bgrat(b0, a0, y0, x0, w1, 15*eps, &ierr1, false);
#ifdef DEBUG_bratio
	REprintf(" ==> new w1=%.15g", *w1);
	if(ierr1) REprintf(" ERROR(code=%d)\n", ierr1) ; else REprintf("\n");
#endif
	if(*w1 == 0 || (0 < *w1 && *w1 < DBL_MIN)) { // w1=0 or very close:
	    // "almost surely" from underflow, try more: [2013-03-04]
// FIXME: it is even better to do this in bgrat *directly* at least for the case
//  !did_bup, i.e., where *w1 = (0 or -Inf) on entry
	    R_ifDEBUG_printf(" denormalized or underflow (?) -> retrying: ");
	    if(did_bup) { // re-do that part on log scale:
		*w1 = bup(b0-n, a0, y0, x0, n, eps, true);
	    }
	    else *w1 = ML_NEGINF; // = 0 on log-scale
	    bgrat(b0, a0, y0, x0, w1, 15*eps, &ierr1, true);
	    if(ierr1) *ierr = 10 + ierr1;
#ifdef DEBUG_bratio
	    REprintf(" ==> new log(w1)=%.15g", *w1);
	    if(ierr1) REprintf(" Error(code=%d)\n", ierr1) ; else REprintf("\n");
#endif
	    goto L_end_from_w1_log;
	}
	// else
	if(ierr1) *ierr = 10 + ierr1;
	if(*w1 < 0)
	    MATHLIB_WARNING4("bratio(a=%g, b=%g, x=%g): bgrat() -> w1 = %g",
			     a,b,x, *w1);
	goto L_end_from_w1;
    }
    else { /* L30: -------------------- both  a, b > 1  {a0 > 1  &  b0 > 1} ---*/

	/* lambda := a y - b x  =  (a + b)y  =  a - (a+b)x    {using x + y == 1},
	 * ------ using the numerically best version : */
	lambda = R_FINITE(a+b)
	    ? ((a > b) ? (a + b) * y - b : a - (a + b) * x)
	    : a*y - b*x;
	do_swap = (lambda < 0.);
	if (do_swap) {
	    lambda = -lambda;
	    SET_0_swap;
	} else {
	    SET_0_noswap;
	}

	R_ifDEBUG_printf("  L30:  both  a, b > 1; |lambda| = %#g, do_swap = %d\n",
			 lambda, do_swap);

	if (b0 < 40.) {
	    R_ifDEBUG_printf("  b0 < 40;");
	    if (b0 * x0 <= 0.7
		|| (log_p && lambda > 650.)) // << added 2010-03; svn r51327
		goto L_w_bpser;
	    else
		goto L140;
	}
	else if (a0 > b0) { /* ----  a0 > b0 >= 40  ---- */
	    R_ifDEBUG_printf("  a0 > b0 >= 40;");
	    if (b0 <= 100. || lambda > b0 * 0.03)
		goto L_bfrac;

	} else if (a0 <= 100.) {
	    R_ifDEBUG_printf("  a0 <= 100; a0 <= b0 >= 40;");
	    goto L_bfrac;
	}
	else if (lambda > a0 * 0.03) {
	    R_ifDEBUG_printf("  b0 >= a0 > 100; lambda > a0 * 0.03 ");
	    goto L_bfrac;
	}

	/* else if none of the above    L180: */
	*w = basym(a0, b0, lambda, eps * 100., log_p);
	*w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
	R_ifDEBUG_printf("  b0 >= a0 > 100; lambda <= a0 * 0.03: *w:= basym(*) =%.15g\n",
			 *w);
	goto L_end;

    } /* else: a, b > 1 */

/*            EVALUATION OF THE APPROPRIATE ALGORITHM */

L_w_bpser: // was L100
    *w = bpser(a0, b0, x0, eps, log_p);
    *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
    R_ifDEBUG_printf(" L_w_bpser: *w := bpser(*) = %.15g\n", *w);
    goto L_end;

L_w1_bpser:  // was L110
    *w1 = bpser(b0, a0, y0, eps, log_p);
    *w  = log_p ? R_Log1_Exp(*w1) : 0.5 - *w1 + 0.5;
    R_ifDEBUG_printf(" L_w1_bpser: *w1 := bpser(*) = %.15g\n", *w1);
    goto L_end;

L_bfrac:
    *w = bfrac(a0, b0, x0, y0, lambda, eps * 15., log_p);
    *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
    R_ifDEBUG_printf(" L_bfrac: *w := bfrac(*) = %g\n", *w);
    goto L_end;

L140:
    /* b0 := fractional_part( b0 )  in (0, 1]  */
    n = (int) b0;
    b0 -= n;
    if (b0 == 0.) {
	--n; b0 = 1.;
    }

    *w = bup(b0, a0, y0, x0, n, eps, false);

    if(*w < DBL_MIN && log_p) { /* do not believe it; try bpser() : */
	R_ifDEBUG_printf(" L140: bup(b0=%g,..)=%.15g < DBL_MIN - not used; ", b0, *w);
	/*revert: */ b0 += n;
	/* which is only valid if b0 <= 1 || b0*x0 <= 0.7 */
	goto L_w_bpser;
    }
    R_ifDEBUG_printf(" L140: *w := bup(b0=%g,..) = %.15g; ", b0, *w);
    if (x0 <= 0.7) {
	/* log_p :  TODO:  w = bup(.) + bpser(.)  -- not so easy to use log-scale */
	*w += bpser(a0, b0, x0, eps, /* log_p = */ false);
	R_ifDEBUG_printf(" x0 <= 0.7: *w := *w + bpser(*) = %.15g\n", *w);
	goto L_end_from_w;
    }
    /* L150: */
    if (a0 <= 15.) {
	n = 20;
	*w += bup(a0, b0, x0, y0, n, eps, false);
	R_ifDEBUG_printf("\n a0 <= 15: *w := *w + bup(*) = %.15g;", *w);
	a0 += n;
    }
    R_ifDEBUG_printf(" bgrat(*, w=%.15g) ", *w);
    bgrat(a0, b0, x0, y0, w, 15*eps, &ierr1, false);
    if(ierr1) *ierr = 10 + ierr1;
#ifdef DEBUG_bratio
    REprintf("==> new w=%.15g", *w);
    if(ierr1) REprintf(" Error(code=%d)\n", ierr1) ; else REprintf("\n");
#endif
    goto L_end_from_w;


/* TERMINATION OF THE PROCEDURE */

L200:
    if (a == 0.) { *ierr = 6;    return; }
    // else:
L201: *w  = R_D__0; *w1 = R_D__1; return;

L210:
    if (b == 0.) { *ierr = 7;    return; }
    // else:
L211: *w  = R_D__1; *w1 = R_D__0; return;

L_end_from_w:
    if(log_p) {
	*w1 = log1p(-*w);
	*w  = log(*w);
    } else {
	*w1 = 0.5 - *w + 0.5;
    }
    goto L_end;

L_end_from_w1:
    if(log_p) {
	*w  = log1p(-*w1);
	*w1 = log(*w1);
    } else {
	*w = 0.5 - *w1 + 0.5;
    }
    goto L_end;

L_end_from_w1_log:
    // *w1 = log(w1) already; w = 1 - w1  ==> log(w) = log(1 - w1) = log(1 - exp(*w1))
    if(log_p) {
	*w = R_Log1_Exp(*w1);
    } else {
	*w  = /* 1 - exp(*w1) */ -expm1(*w1);
	*w1 = exp(*w1);
    }
    goto L_end;


L_end:
    if (do_swap) { /* swap */
	double t = *w; *w = *w1; *w1 = t;
    }
    return;

} /* bratio */




//pbeta_raw
double pbeta_raw(double x, double a, double b, int lower_tail, int log_p)
{
    // treat limit cases correctly here:
    if(a == 0 || b == 0 || !R_FINITE(a) || !R_FINITE(b)) {
	// NB:  0 < x < 1 :
	if(a == 0 && b == 0) // point mass 1/2 at each of {0,1} :
	    return (log_p ? -M_LN2 : 0.5);
	if (a == 0 || a/b == 0) // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
	    return R_DT_1;
	if (b == 0 || b/a == 0) // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
	    return R_DT_0;
	// else, remaining case:  a = b = Inf : point mass 1 at 1/2
	if (x < 0.5) return R_DT_0; else return R_DT_1;
    }
    // Now:  0 < a < Inf;  0 < b < Inf

    double x1 = 0.5 - x + 0.5, w, wc;
    int ierr;
    //====
    bratio(a, b, x, x1, &w, &wc, &ierr, log_p); /* -> ./toms708.c */
    //====
    // ierr in {10,14} <==> bgrat() error code ierr-10 in 1:4; for 1 and 4, warned *there*
    if(ierr && ierr != 11 && ierr != 14)
	MATHLIB_WARNING4("pbeta_raw(%g, a=%g, b=%g, ..) -> bratio() gave error code %d", x, a,b, ierr);
    return lower_tail ? w : wc;
} /* pbeta_raw() */



//pbeta
double pbeta(double x, double a, double b, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a) || ISNAN(b)) return x + a + b;
#endif

    if (a < 0 || b < 0) ML_ERR_return_NAN;
    // allowing a==0 and b==0  <==> treat as one- or two-point mass

    if (x <= 0)
	return R_DT_0;
    if (x >= 1)
	return R_DT_1;

    return pbeta_raw(x, a, b, lower_tail, log_p);
}


//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
//pt
//////////////////////////////////////////////////////////////////////////////

//pt from pt.c
double pt(double x, double n, int lower_tail, int log_p)
{
/* return  P[ T <= x ]	where
 * T ~ t_{n}  (t distrib. with n degrees of freedom).

 *	--> ./pnt.c for NON-central
 */
    double val, nx;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
	return x + n;
#endif
    if (n <= 0.0) ML_ERR_return_NAN;

    if(!R_FINITE(x))
	return (x < 0) ? R_DT_0 : R_DT_1;
    if(!R_FINITE(n))
	return pnorm(x, 0.0, 1.0, lower_tail, log_p);

#ifdef R_version_le_260
    if (n > 4e5) { /*-- Fixme(?): test should depend on `n' AND `x' ! */
	/* Approx. from	 Abramowitz & Stegun 26.7.8 (p.949) */
	val = 1./(4.*n);
	return pnorm(x*(1. - val)/sqrt(1. + x*x*2.*val), 0.0, 1.0,
		     lower_tail, log_p);
    }
#endif

    nx = 1 + (x/n)*x;
    /* FIXME: This test is probably losing rather than gaining precision,
     * now that pbeta(*, log_p = TRUE) is much better.
     * Note however that a version of this test *is* needed for x*x > D_MAX */
    if(nx > 1e100) { /* <==>  x*x > 1e100 * n  */
	/* Danger of underflow. So use Abramowitz & Stegun 26.5.4
	   pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
	   with z = 1/nx,  a = n/2,  b= 1/2 :
	*/
	double lval;
	lval = -0.5*n*(2*log(fabs(x)) - log(n))
		- lbeta(0.5*n, 0.5) - log(0.5*n);
	val = log_p ? lval : exp(lval);
    } else {
	val = (n > x * x)
	    ? pbeta (x * x / (n + x * x), 0.5, n / 2., /*lower_tail*/0, log_p)
	    : pbeta (1. / nx,             n / 2., 0.5, /*lower_tail*/1, log_p);
    }

    /* Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):*/
    if(x <= 0.)
	lower_tail = !lower_tail;

    if(log_p) {
	if(lower_tail) return log1p(-0.5*exp(val));
	else return val - M_LN2; /* = log(.5* pbeta(....)) */
    }
    else {
	val /= 2.;
	return R_D_Cval(val);
    }
}



//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// exposed functions
//////////////////////////////////////////////////////////////////////////////

//cdf of t dist, no log
double cdf_t(double x, double n, int lower_tail){
    return pt(x, n, lower_tail, false);
}

//cdf of t, with log p
double cdf_t_log(double x, double n, int lower_tail){
    return pt(x, n, lower_tail, true);
}

//chisq cdf, no log 
double cdf_chisq(double x, double df, int lower_tail){
    return pchisq(x, df, lower_tail, false);
}

//chisq cdf, with log p
double cdf_chisq_log(double x, double df, int lower_tail){
    return pchisq(x, df, lower_tail, true);
}

//////////////////////////////////////////////////////////////////////////////
