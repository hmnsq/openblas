#अन्तर्भाव <math.h>
#अन्तर्भाव <float.h>
#अन्तर्भाव "common.h"
#यदि_परिभाषितम् FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS

शून्य NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#अन्यथा

शून्य CNAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#यदि_अन्तः

#यदि_परिभाषितम् DOUBLE
  दीर्घ द्विगुण safmin = DBL_MIN;
#अन्यथा
  दीर्घ द्विगुण safmin = FLT_MIN;
#यदि_अन्तः

#यदि defined(__i386__) || defined(__x86_64__) || defined(__ia64__) || defined(_M_X64) || defined(_M_IX86)

  दीर्घ द्विगुण da = *DA;
  दीर्घ द्विगुण db = *DB;
  दीर्घ द्विगुण c;
  दीर्घ द्विगुण s;
  दीर्घ द्विगुण r, z;
  दीर्घ द्विगुण sigma, dascal,dbscal;

  दीर्घ द्विगुण ada = fabsl(da);
  दीर्घ द्विगुण adb = fabsl(db);
  दीर्घ द्विगुण maxab = MAX(ada,adb);
  दीर्घ द्विगुण safmax;
  दीर्घ द्विगुण scale;


#यदि_अपरिभाषितम् CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#यदि_अन्तः

  यदि (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } अन्यथा यदि (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } अन्यथा {
  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);
    यदि (ada > adb)
	sigma = copysign(1.,da);
    अन्यथा
	sigma = copysign(1.,db);
    dascal = da / scale;
    dbscal = db / scale;
    r = sigma * (scale * वर्गमूलम्(dascal * dascal + dbscal * dbscal));
    c = da / r;
    s = db / r;
    z = ONE;
    यदि (ada > adb) z = s;
    यदि ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }

#अन्यथा
  FLOAT da = *DA;
  FLOAT db = *DB;
  FLOAT c  = *C;
  FLOAT s  = *S;
  FLOAT sigma;
  FLOAT r, z;

  FLOAT ada = fabs(da);
  FLOAT adb = fabs(db);
  FLOAT maxab = MAX(ada,adb);
  दीर्घ द्विगुण safmax ;
  FLOAT scale ;

  safmax = 1./safmin;
  scale = MIN(MAX(safmin,maxab), safmax);

  यदि (ada > adb)
	sigma = copysign(1.,da);
    अन्यथा
	sigma = copysign(1.,db);

#यदि_अपरिभाषितम् CBLAS
  PRINT_DEBUG_NAME;
#अन्यथा
  PRINT_DEBUG_CNAME;
#यदि_अन्तः


  यदि (adb == ZERO) {
    *C = ONE;
    *S = ZERO;
    *DB = ZERO;
  } अन्यथा यदि (ada == ZERO) {
    *C = ZERO;
    *S = ONE;
    *DA = *DB;
    *DB = ONE;
  } अन्यथा {
    FLOAT aa = da / scale;
    FLOAT bb = db / scale;

    r = sigma * scale * वर्गमूलम्(aa * aa + bb * bb);
    c = da / r;
    s = db / r;
    z = ONE;
    यदि (ada > adb) z = s;
    यदि ((ada <= adb) && (c != ZERO)) z = ONE / c;

    *C = c;
    *S = s;
    *DA = r;
    *DB = z;
  }
#यदि_अन्तः

  प्रत्यागम;
}
