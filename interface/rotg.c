#ಅನ್ತರ್ಭಾವ <math.h>
#ಅನ್ತರ್ಭಾವ <float.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long ದ್ವಿಗುಣ safmin = DBL_MIN;
#ಅನ್ಯಥಾ
  long ದ್ವಿಗುಣ safmin = FLT_MIN;
#endif

#ಯದಿ defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long ದ್ವಿಗುಣ da = *DA;
  long ದ್ವಿಗುಣ db = *DB;
  long ದ್ವಿಗುಣ c;
  long ದ್ವಿಗುಣ s;
  long ದ್ವಿಗುಣ r, z;
  long ದ್ವಿಗುಣ sigma, dascal,dbscal;

  long ದ್ವಿಗುಣ ada = fabsl(da);
  long ದ್ವಿಗುಣ adb = fabsl(db);
  long ದ್ವಿಗುಣ maxab = MAX(ada,adb);
  long ದ್ವಿಗುಣ safmax;
  long ದ್ವಿಗುಣ scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ಅನ್ಯಥಾ
  PRINT_DEBUG_CNAME;
#endif

  ಯದಿ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ಅನ್ಯಥಾ ಯದಿ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ಅನ್ಯಥಾ {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    ಯದಿ (ada > adb)
	sigma = copysign(1.,da);
    ಅನ್ಯಥಾ
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    ಯದಿ (ada > adb) z = s;
    ಯದಿ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#ಅನ್ಯಥಾ
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long ದ್ವಿಗುಣ safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  ಯದಿ (ada > adb)
	sigma = copysign(1.,da);
    ಅನ್ಯಥಾ
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ಅನ್ಯಥಾ
  PRINT_DEBUG_CNAME;
#endif


  ಯದಿ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ಅನ್ಯಥಾ ಯದಿ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ಅನ್ಯಥಾ {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    ಯದಿ (ada > adb) z = s;
    ಯದಿ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  ಪ್ರತ್ಯಾಗಮ;
}
