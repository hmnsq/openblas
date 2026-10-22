#شامل <math.h>
#شامل <float.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ورنہ

خالی CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#اگر_ختم

#اگر_تعریف_شدہ DOUBLE
  طویل دوہرا safmin = DBL_MIN;
#ورنہ
  طویل دوہرا safmin = FLT_MIN;
#اگر_ختم

#اگر defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  طویل دوہرا da = *DA;
  طویل دوہرا db = *DB;
  طویل دوہرا c;
  طویل دوہرا s;
  طویل دوہرا r, z;
  طویل دوہرا sigma, dascal,dbscal;

  طویل دوہرا ada = fabsl(da);
  طویل دوہرا adb = fabsl(db);
  طویل دوہرا maxab = MAX(ada,adb);
  طویل دوہرا safmax;
  طویل دوہرا scale;


#اگر_غیر_تعریف_شدہ CBLAS
  PRINT_DEBUG_NAME;
#ورنہ
  PRINT_DEBUG_CNAME;
#اگر_ختم

  اگر (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ورنہ اگر (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ورنہ {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    اگر (ada > adb)
	sigma = copysign(1.,da);
    ورنہ
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * جذر(dascal * dascal + dbscal * dbscal));
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

#ورنہ
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  طویل دوہرا safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  اگر (ada > adb)
	sigma = copysign(1.,da);
    ورنہ
	sigma = copysign(1.,db);

#اگر_غیر_تعریف_شدہ CBLAS
  PRINT_DEBUG_NAME;
#ورنہ
  PRINT_DEBUG_CNAME;
#اگر_ختم


  اگر (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ورنہ اگر (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ورنہ {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * جذر(aa * aa + bb * bb);
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
#اگر_ختم

  واپس;
}
