#ضمن <math.h>
#ضمن <float.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#אחרת

فراغ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long مزدوج safmin = DBL_MIN;
#אחרת
  long مزدوج safmin = FLT_MIN;
#endif

#אם defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long مزدوج da = *DA;
  long مزدوج db = *DB;
  long مزدوج c;
  long مزدوج s;
  long مزدوج r, z;
  long مزدوج sigma, dascal,dbscal;

  long مزدوج ada = fabsl(da);
  long مزدوج adb = fabsl(db);
  long مزدوج maxab = MAX(ada,adb);
  long مزدوج safmax;
  long مزدوج scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#אחרת
  PRINT_DEBUG_CNAME;
#endif

  אם (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } אחרת אם (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } אחרת {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    אם (ada > adb)
	sigma = copysign(1.,da);
    אחרת
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    אם (ada > adb) z = s;
    אם ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#אחרת
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long مزدوج safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  אם (ada > adb)
	sigma = copysign(1.,da);
    אחרת
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#אחרת
  PRINT_DEBUG_CNAME;
#endif


  אם (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } אחרת אם (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } אחרת {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    אם (ada > adb) z = s;
    אם ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  החזר;
}
