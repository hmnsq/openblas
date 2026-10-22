#അന്തര്ഭാവ <math.h>
#അന്തര്ഭാവ <float.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ശൂന്യം NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#അല്ലെങ്കിൽ

ശൂന്യം CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long ഇരട്ട safmin = DBL_MIN;
#അല്ലെങ്കിൽ
  long ഇരട്ട safmin = FLT_MIN;
#endif

#എങ്കിൽ defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long ഇരട്ട da = *DA;
  long ഇരട്ട db = *DB;
  long ഇരട്ട c;
  long ഇരട്ട s;
  long ഇരട്ട r, z;
  long ഇരട്ട sigma, dascal,dbscal;

  long ഇരട്ട ada = fabsl(da);
  long ഇരട്ട adb = fabsl(db);
  long ഇരട്ട maxab = MAX(ada,adb);
  long ഇരട്ട safmax;
  long ഇരട്ട scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#അല്ലെങ്കിൽ
  PRINT_DEBUG_CNAME;
#endif

  എങ്കിൽ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } അല്ലെങ്കിൽ എങ്കിൽ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } അല്ലെങ്കിൽ {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    എങ്കിൽ (ada > adb)
	sigma = copysign(1.,da);
    അല്ലെങ്കിൽ
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    എങ്കിൽ (ada > adb) z = s;
    എങ്കിൽ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#അല്ലെങ്കിൽ
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long ഇരട്ട safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  എങ്കിൽ (ada > adb)
	sigma = copysign(1.,da);
    അല്ലെങ്കിൽ
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#അല്ലെങ്കിൽ
  PRINT_DEBUG_CNAME;
#endif


  എങ്കിൽ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } അല്ലെങ്കിൽ എങ്കിൽ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } അല്ലെങ്കിൽ {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    എങ്കിൽ (ada > adb) z = s;
    എങ്കിൽ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  തിരികെ;
}
