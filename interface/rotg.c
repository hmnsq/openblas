#شامل <math.h>
#شامل <float.h>
#شامل "common.h"
#ifdef FUNCTION_PROFILE
#شامل "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

تهی NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#وگرنه

تهی CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long دوگانه safmin = DBL_MIN;
#وگرنه
  long دوگانه safmin = FLT_MIN;
#endif

#اگر defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long دوگانه da = *DA;
  long دوگانه db = *DB;
  long دوگانه c;
  long دوگانه s;
  long دوگانه r, z;
  long دوگانه sigma, dascal,dbscal;

  long دوگانه ada = fabsl(da);
  long دوگانه adb = fabsl(db);
  long دوگانه maxab = MAX(ada,adb);
  long دوگانه safmax;
  long دوگانه scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#وگرنه
  PRINT_DEBUG_CNAME;
#endif

  اگر (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } وگرنه اگر (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } وگرنه {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    اگر (ada > adb)
	sigma = copysign(1.,da);
    وگرنه
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    اگر (ada > adb) z = s;
    اگر ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#وگرنه
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long دوگانه safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  اگر (ada > adb)
	sigma = copysign(1.,da);
    وگرنه
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#وگرنه
  PRINT_DEBUG_CNAME;
#endif


  اگر (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } وگرنه اگر (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } وگرنه {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    اگر (ada > adb) z = s;
    اگر ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  بازگشت;
}
