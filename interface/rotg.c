#ଅନ୍ତର୍ଭାଵ <math.h>
#ଅନ୍ତର୍ଭାଵ <float.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ଶୂନ୍ୟ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long ଦ୍ୱିଗୁଣ safmin = DBL_MIN;
#ଅନ୍ୟଥା
  long ଦ୍ୱିଗୁଣ safmin = FLT_MIN;
#endif

#ଯଦି defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long ଦ୍ୱିଗୁଣ da = *DA;
  long ଦ୍ୱିଗୁଣ db = *DB;
  long ଦ୍ୱିଗୁଣ c;
  long ଦ୍ୱିଗୁଣ s;
  long ଦ୍ୱିଗୁଣ r, z;
  long ଦ୍ୱିଗୁଣ sigma, dascal,dbscal;

  long ଦ୍ୱିଗୁଣ ada = fabsl(da);
  long ଦ୍ୱିଗୁଣ adb = fabsl(db);
  long ଦ୍ୱିଗୁଣ maxab = MAX(ada,adb);
  long ଦ୍ୱିଗୁଣ safmax;
  long ଦ୍ୱିଗୁଣ scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ଅନ୍ୟଥା
  PRINT_DEBUG_CNAME;
#endif

  ଯଦି (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ଅନ୍ୟଥା ଯଦି (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ଅନ୍ୟଥା {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    ଯଦି (ada > adb)
	sigma = copysign(1.,da);
    ଅନ୍ୟଥା
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    ଯଦି (ada > adb) z = s;
    ଯଦି ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#ଅନ୍ୟଥା
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long ଦ୍ୱିଗୁଣ safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  ଯଦି (ada > adb)
	sigma = copysign(1.,da);
    ଅନ୍ୟଥା
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ଅନ୍ୟଥା
  PRINT_DEBUG_CNAME;
#endif


  ଯଦି (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ଅନ୍ୟଥା ଯଦି (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ଅନ୍ୟଥା {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    ଯଦି (ada > adb) z = s;
    ଯଦି ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  ଫେରନ୍ତୁ;
}
