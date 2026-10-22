#ضمن <math.h>
#ضمن <float.h>
#ضمن "common.h"
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#وإلا

فراغ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#جے_انت

#جے_تعریف_شدہ DOUBLE
  لما مزدوج safmin = DBL_MIN;
#وإلا
  لما مزدوج safmin = FLT_MIN;
#جے_انت

#إذا defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  لما مزدوج da = *DA;
  لما مزدوج db = *DB;
  لما مزدوج c;
  لما مزدوج s;
  لما مزدوج r, z;
  لما مزدوج sigma, dascal,dbscal;

  لما مزدوج ada = fabsl(da);
  لما مزدوج adb = fabsl(db);
  لما مزدوج maxab = MAX(ada,adb);
  لما مزدوج safmax;
  لما مزدوج scale;


#جے_بنا_تعریف CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#جے_انت

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
  لما مزدوج safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  إذا (ada > adb)
	sigma = copysign(1.,da);
    وإلا
	sigma = copysign(1.,db);

#جے_بنا_تعریف CBLAS
  PRINT_DEBUG_NAME;
#وإلا
  PRINT_DEBUG_CNAME;
#جے_انت


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
#جے_انت

  أرجع;
}
