#अन्तर्भाव <math.h>
#अन्तर्भाव <float.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif

OPENBLAS_EXPORT

#ifndef CBLAS

शून्य NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#अञ्ञथा

शून्य CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#endif

#ifdef DOUBLE
  long द्विगुण safmin = DBL_MIN;
#अञ्ञथा
  long द्विगुण safmin = FLT_MIN;
#endif

#सचे defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  long द्विगुण da = *DA;
  long द्विगुण db = *DB;
  long द्विगुण c;
  long द्विगुण s;
  long द्विगुण r, z;
  long द्विगुण sigma, dascal,dbscal;

  long द्विगुण ada = fabsl(da);
  long द्विगुण adb = fabsl(db);
  long द्विगुण maxab = MAX(ada,adb);
  long द्विगुण safmax;
  long द्विगुण scale;


#ifndef CBLAS
  PRINT_DEBUG_NAME;
#अञ्ञथा
  PRINT_DEBUG_CNAME;
#endif

  सचे (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } अञ्ञथा सचे (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } अञ्ञथा {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    सचे (ada > adb)
	sigma = copysign(1.,da);
    अञ्ञथा
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * sqrt(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    सचे (ada > adb) z = s;
    सचे ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#अञ्ञथा
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  long द्विगुण safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  सचे (ada > adb)
	sigma = copysign(1.,da);
    अञ्ञथा
	sigma = copysign(1.,db);

#ifndef CBLAS
  PRINT_DEBUG_NAME;
#अञ्ञथा
  PRINT_DEBUG_CNAME;
#endif


  सचे (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } अञ्ञथा सचे (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } अञ्ञथा {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * sqrt(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    सचे (ada > adb) z = s;
    सचे ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#endif

  पटिनिवत्त;
}
