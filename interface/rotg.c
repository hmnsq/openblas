#ಅನ್ತರ್ಭಾವ <math.h>
#ಅನ್ತರ್ಭಾವ <float.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS

ಶೂನ್ಯ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DOUBLE
  ದೀರ್ಘ ದ್ವಿಗುಣ safmin = DBL_MIN;
#ಇಲ್ಲದಿದ್ದರೆ
  ದೀರ್ಘ ದ್ವಿಗುಣ safmin = FLT_MIN;
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  ದೀರ್ಘ ದ್ವಿಗುಣ da = *DA;
  ದೀರ್ಘ ದ್ವಿಗುಣ db = *DB;
  ದೀರ್ಘ ದ್ವಿಗುಣ c;
  ದೀರ್ಘ ದ್ವಿಗುಣ s;
  ದೀರ್ಘ ದ್ವಿಗುಣ r, z;
  ದೀರ್ಘ ದ್ವಿಗುಣ sigma, dascal,dbscal;

  ದೀರ್ಘ ದ್ವಿಗುಣ ada = fabsl(da);
  ದೀರ್ಘ ದ್ವಿಗುಣ adb = fabsl(db);
  ದೀರ್ಘ ದ್ವಿಗುಣ maxab = MAX(ada,adb);
  ದೀರ್ಘ ದ್ವಿಗುಣ safmax;
  ದೀರ್ಘ ದ್ವಿಗುಣ scale;


#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
  PRINT_DEBUG_NAME;
#ಇಲ್ಲದಿದ್ದರೆ
  PRINT_DEBUG_CNAME;
#ಆದರೆ_ಅಂತ್ಯ

  ಆದರೆ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ಇಲ್ಲದಿದ್ದರೆ {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    ಆದರೆ (ada > adb)
	sigma = copysign(1.,da);
    ಇಲ್ಲದಿದ್ದರೆ
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * ವರ್ಗಮೂಲ(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    ಆದರೆ (ada > adb) z = s;
    ಆದರೆ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#ಇಲ್ಲದಿದ್ದರೆ
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  ದೀರ್ಘ ದ್ವಿಗುಣ safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  ಆದರೆ (ada > adb)
	sigma = copysign(1.,da);
    ಇಲ್ಲದಿದ್ದರೆ
	sigma = copysign(1.,db);

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
  PRINT_DEBUG_NAME;
#ಇಲ್ಲದಿದ್ದರೆ
  PRINT_DEBUG_CNAME;
#ಆದರೆ_ಅಂತ್ಯ


  ಆದರೆ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ಇಲ್ಲದಿದ್ದರೆ {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * ವರ್ಗಮೂಲ(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    ಆದರೆ (ada > adb) z = s;
    ಆದರೆ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#ಆದರೆ_ಅಂತ್ಯ

  ಹಿಂತಿರುಗಿ;
}
