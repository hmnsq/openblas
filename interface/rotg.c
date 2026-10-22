#ضمن <math.h>
#ضمن <float.h>
#ضمن "common.h"
#اگر_تعریف_شده FUNCTION_PROFILE
#ضمن "functable.h"
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS

فراغ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#وإلا

فراغ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#پایان_اگر

#اگر_تعریف_شده DOUBLE
  دراز مزدوج safmin = DBL_MIN;
#وإلا
  دراز مزدوج safmin = FLT_MIN;
#پایان_اگر

#إذا defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  دراز مزدوج da = *DA;
  دراز مزدوج db = *DB;
  دراز مزدوج c;
  دراز مزدوج s;
  دراز مزدوج r, z;
  دراز مزدوج sigma, dascal,dbscal;

  دراز مزدوج ada = fabsl(da);
  دراز مزدوج adb = fabsl(db);
  دراز مزدوج maxab = MAX(ada,adb);
  دراز مزدوج safmax;
  دراز مزدوج scale;


#اگر_تعریف_نشده CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#پایان_اگر

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
    r = sigma * (scale * جذر(dascal * dascal + dbscal * dbscal));
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
  دراز مزدوج safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  إذا (ada > adb)
	sigma = copysign(1.,da);
    وإلا
	sigma = copysign(1.,db);

#اگر_تعریف_نشده CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#پایان_اگر


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

    r = sigma * scale * جذر(aa * aa + bb * bb);
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
#پایان_اگر

  أرجع;
}
