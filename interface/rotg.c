#અન્તર્ભાવ <math.h>
#અન્તર્ભાવ <float.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#નહીંતર

શૂન્ય CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long દ્વિગુણ safmin = DBL_MIN;
#નહીંતર
  long દ્વિગુણ safmin = FLT_MIN;
#endif

#જો defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long દ્વિગુણ da = *DA;
  long દ્વિગુણ db = *DB;
  long દ્વિગુણ c;
  long દ્વિગુણ s;
  long દ્વિગુણ r, z;
  long દ્વિગુણ sigma, dascal,dbscal;

  long દ્વિગુણ ada = fabsl(da);
  long દ્વિગુણ adb = fabsl(db);
  long દ્વિગુણ maxab = MAX(ada,adb);
  long દ્વિગુણ safmax;
  long દ્વિગુણ scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#નહીંતર
  PRINT_DEBUG_CNAME;
#endif

  જો (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } નહીંતર જો (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } નહીંતર {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    જો (ada > adb)
	sigma = copysign(1.,da);
    નહીંતર
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    જો (ada > adb) z = s;
    જો ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#નહીંતર
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long દ્વિગુણ safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  જો (ada > adb)
	sigma = copysign(1.,da);
    નહીંતર
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#નહીંતર
  PRINT_DEBUG_CNAME;
#endif


  જો (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } નહીંતર જો (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } નહીંતર {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    જો (ada > adb) z = s;
    જો ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  પાછા;
}
