#include <math.h>
#include <float.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

void NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ບໍ່ດັ່ງນັ້ນ

void CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long dviguu_na safmin = DBL_MIN;
#ບໍ່ດັ່ງນັ້ນ
  long dviguu_na safmin = FLT_MIN;
#endif

#ຖ້າ defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long dviguu_na da = *DA;
  long dviguu_na db = *DB;
  long dviguu_na c;
  long dviguu_na s;
  long dviguu_na r, z;
  long dviguu_na sigma, dascal,dbscal;

  long dviguu_na ada = fabsl(da);
  long dviguu_na adb = fabsl(db);
  long dviguu_na maxab = MAX(ada,adb);
  long dviguu_na safmax;
  long dviguu_na scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ບໍ່ດັ່ງນັ້ນ
  PRINT_DEBUG_CNAME;
#endif

  ຖ້າ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ບໍ່ດັ່ງນັ້ນ ຖ້າ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ບໍ່ດັ່ງນັ້ນ {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    ຖ້າ (ada > adb)
	sigma = copysign(1.,da);
    ບໍ່ດັ່ງນັ້ນ
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    ຖ້າ (ada > adb) z = s;
    ຖ້າ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#ບໍ່ດັ່ງນັ້ນ
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long dviguu_na safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  ຖ້າ (ada > adb)
	sigma = copysign(1.,da);
    ບໍ່ດັ່ງນັ້ນ
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#ບໍ່ດັ່ງນັ້ນ
  PRINT_DEBUG_CNAME;
#endif


  ຖ້າ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } ບໍ່ດັ່ງນັ້ນ ຖ້າ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } ບໍ່ດັ່ງນັ້ນ {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    ຖ້າ (ada > adb) z = s;
    ຖ້າ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  ສົ່ງຄືນ;
}
