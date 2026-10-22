#అన్తర్భావ <math.h>
#అన్తర్భావ <float.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS

శూన్య NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#లేకపోతే

శూన్య CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#అయితే_ముగింపు

#నిర్వచితమైతే DOUBLE
  దీర్ఘ ద్విగుణ safmin = DBL_MIN;
#లేకపోతే
  దీర్ఘ ద్విగుణ safmin = FLT_MIN;
#అయితే_ముగింపు

#ఒకవేళ defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  దీర్ఘ ద్విగుణ da = *DA;
  దీర్ఘ ద్విగుణ db = *DB;
  దీర్ఘ ద్విగుణ c;
  దీర్ఘ ద్విగుణ s;
  దీర్ఘ ద్విగుణ r, z;
  దీర్ఘ ద్విగుణ sigma, dascal,dbscal;

  దీర్ఘ ద్విగుణ ada = fabsl(da);
  దీర్ఘ ద్విగుణ adb = fabsl(db);
  దీర్ఘ ద్విగుణ maxab = MAX(ada,adb);
  దీర్ఘ ద్విగుణ safmax;
  దీర్ఘ ద్విగుణ scale;


#నిర్వచితం_కాకపోతే CBLAS
  PRINT_DEBUG_NAME;
#లేకపోతే
  PRINT_DEBUG_CNAME;
#అయితే_ముగింపు

  ఒకవేళ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } లేకపోతే ఒకవేళ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } లేకపోతే {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    ఒకవేళ (ada > adb)
	sigma = copysign(1.,da);
    లేకపోతే
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * వర్గమూలం(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    ఒకవేళ (ada > adb) z = s;
    ఒకవేళ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#లేకపోతే
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  దీర్ఘ ద్విగుణ safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  ఒకవేళ (ada > adb)
	sigma = copysign(1.,da);
    లేకపోతే
	sigma = copysign(1.,db);

#నిర్వచితం_కాకపోతే CBLAS
  PRINT_DEBUG_NAME;
#లేకపోతే
  PRINT_DEBUG_CNAME;
#అయితే_ముగింపు


  ఒకవేళ (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } లేకపోతే ఒకవేళ (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } లేకపోతే {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * వర్గమూలం(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    ఒకవేళ (ada > adb) z = s;
    ఒకవేళ ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#అయితే_ముగింపు

  తిరిగి;
}
