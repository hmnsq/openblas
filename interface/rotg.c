#ضمن <math.h>
#ضمن <float.h>
#ضمن "common.h"
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#وإلا

فراغ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long مزدوج safmin = DBL_MIN;
#وإلا
  long مزدوج safmin = FLT_MIN;
#endif

#إذا defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

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
  اطبع_DEBUG_NAME;
#وإلا
  اطبع_DEBUG_CNAME;
#endif

  إذا (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } وإلا إذا (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } وإلا {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    إذا (ada > adb)
	sigma = copysign(1.,da);
    وإلا
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    إذا (ada > adb) z = s;
    إذا ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#وإلا
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

  إذا (ada > adb)
	sigma = copysign(1.,da);
    وإلا
	sigma = copysign(1.,db);

#ifndef CBLAS
  اطبع_DEBUG_NAME;
#وإلا
  اطبع_DEBUG_CNAME;
#endif


  إذا (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } وإلا إذا (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } وإلا {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    إذا (ada > adb) z = s;
    إذا ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  أرجع;
}
