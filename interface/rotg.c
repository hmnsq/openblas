#include <math.h>
#include <float.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#инако

void CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long double safmin = DBL_MIN;
#инако
  long double safmin = FLT_MIN;
#endif

#ако defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long double da = *DA;
  long double db = *DB;
  long double c;
  long double s;
  long double r, z;
  long double sigma, dascal,dbscal;

  long double ada = fabsl(da);
  long double adb = fabsl(db);
  long double maxab = MAX(ada,adb);
  long double safmax;
  long double scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#инако
  PRINT_DEBUG_CNAME;
#endif

  ако (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } инако ако (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } инако {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    ако (ada > adb)
	sigma = copysign(1.,da);
    инако
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    ако (ada > adb) z = s;
    ако ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#инако
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long double safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  ако (ada > adb)
	sigma = copysign(1.,da);
    инако
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#инако
  PRINT_DEBUG_CNAME;
#endif


  ако (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } инако ако (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } инако {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    ако (ada > adb) z = s;
    ако ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  врати;
}
