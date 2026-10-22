#ਅਨ੍ਤਰ੍ਭਾਵ <math.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <float.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
  ਲੰਮਾ ਦੁੱਗਣਾ safmin = DBL_MIN;
#ਨਹੀਂਤਾਂ
  ਲੰਮਾ ਦੁੱਗਣਾ safmin = FLT_MIN;
#ਜੇ_ਅੰਤ

#ਜੇ defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  ਲੰਮਾ ਦੁੱਗਣਾ da = *DA;
  ਲੰਮਾ ਦੁੱਗਣਾ db = *DB;
  ਲੰਮਾ ਦੁੱਗਣਾ c;
  ਲੰਮਾ ਦੁੱਗਣਾ s;
  ਲੰਮਾ ਦੁੱਗਣਾ r, z;
  ਲੰਮਾ ਦੁੱਗਣਾ sigma, dascal,dbscal;

  ਲੰਮਾ ਦੁੱਗਣਾ ada = fabsl(da);
  ਲੰਮਾ ਦੁੱਗਣਾ adb = fabsl(db);
  ਲੰਮਾ ਦੁੱਗਣਾ maxab = MAX(ada,adb);
  ਲੰਮਾ ਦੁੱਗਣਾ safmax;
  ਲੰਮਾ ਦੁੱਗਣਾ scale;


#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS
  PRINT_DEBUG_NAME;
#ਨਹੀਂਤਾਂ
  PRINT_DEBUG_CNAME;
#ਜੇ_ਅੰਤ

  ਜੇ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ਨਹੀਂਤਾਂ ਜੇ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ਨਹੀਂਤਾਂ {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    ਜੇ (ada > adb)
	sigma = copysign(1.,da);
    ਨਹੀਂਤਾਂ
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * ਵਰਗਮੂਲ(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    ਜੇ (ada > adb) z = s;
    ਜੇ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#ਨਹੀਂਤਾਂ
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  ਲੰਮਾ ਦੁੱਗਣਾ safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  ਜੇ (ada > adb)
	sigma = copysign(1.,da);
    ਨਹੀਂਤਾਂ
	sigma = copysign(1.,db);

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS
  PRINT_DEBUG_NAME;
#ਨਹੀਂਤਾਂ
  PRINT_DEBUG_CNAME;
#ਜੇ_ਅੰਤ


  ਜੇ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ਨਹੀਂਤਾਂ ਜੇ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ਨਹੀਂਤਾਂ {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * ਵਰਗਮੂਲ(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    ਜੇ (ada > adb) z = s;
    ਜੇ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#ਜੇ_ਅੰਤ

  ਵਾਪਸ;
}
