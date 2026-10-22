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
ਖਾਲੀ CNAME(ਖਾਲੀ *VDA, ਖਾਲੀ *VDB, FLOAT *C, ਖਾਲੀ *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#ਜੇ_ਅੰਤ /* CBLAS */

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DOUBLE
  ਲੰਮਾ ਦੁੱਗਣਾ safmin = DBL_MIN;
  ਲੰਮਾ ਦੁੱਗਣਾ rtmin = ਵਰਗਮੂਲ(DBL_MIN/DBL_EPSILON);
#ਨਹੀਂਤਾਂ
  ਲੰਮਾ ਦੁੱਗਣਾ safmin = FLT_MIN;
  ਲੰਮਾ ਦੁੱਗਣਾ rtmin = ਵਰਗਮੂਲ(FLT_MIN/FLT_EPSILON);
#ਜੇ_ਅੰਤ


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  ਲੰਮਾ ਦੁੱਗਣਾ d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();
  ਲੰਮਾ ਦੁੱਗਣਾ safmax = 1./safmin;
  safmax *= 2.;
  ਕਰੋ {
	ਜੇ (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		ਰੋਕੋ;
	}

#ਜੇ defined DOUBLE
	ਲੰਮਾ ਦੁੱਗਣਾ rtmax = safmax /DBL_EPSILON;
#ਨਹੀਂਤਾਂ
	ਲੰਮਾ ਦੁੱਗਣਾ rtmax = safmax /FLT_EPSILON;
#ਜੇ_ਅੰਤ

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	ਜੇ (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		ਜੇ (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ਰੋਕੋ;
		} ਨਹੀਂਤਾਂ ਜੇ ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ਰੋਕੋ;
		} ਨਹੀਂਤਾਂ {
			ਲੰਮਾ ਦੁੱਗਣਾ g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =ਵਰਗਮੂਲ(safmax/2.);
			ਜੇ (g1 > rtmin && g1 < rtmax) { // unscaled
				d = ਵਰਗਮੂਲ(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				ਰੋਕੋ;
			} ਨਹੀਂਤਾਂ { // scaled algorithm
				ਲੰਮਾ ਦੁੱਗਣਾ u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = ਵਰਗਮੂਲ ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				ਰੋਕੋ;
			}
		}

	} ਨਹੀਂਤਾਂ {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = ਵਰਗਮੂਲ(safmax / 4.);
		ਜੇ ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			ਲੰਮਾ ਦੁੱਗਣਾ h = ada + adb;
			ਦੁੱਗਣਾ adahsq = ਵਰਗਮੂਲ(ada * h);
			ਜੇ (ada >= h *safmin) {
				*C = ਵਰਗਮੂਲ(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				ਜੇ ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} ਨਹੀਂਤਾਂ {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} ਨਹੀਂਤਾਂ {
				*C = ada / adahsq;
				ਜੇ (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} ਨਹੀਂਤਾਂ {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			ਰੋਕੋ;
		} ਨਹੀਂਤਾਂ {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			ਲੰਮਾ ਦੁੱਗਣਾ v,w,f2,g2,h;
			ਲੰਮਾ ਦੁੱਗਣਾ u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			ਜੇ (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} ਨਹੀਂਤਾਂ { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			ਜੇ ( f2 >= h * safmin) {
				*C = ਵਰਗਮੂਲ ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				ਜੇ ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=ਵਰਗਮੂਲ(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} ਨਹੀਂਤਾਂ {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} ਨਹੀਂਤਾਂ { // intermediates might overflow
				d = ਵਰਗਮੂਲ ( f2 * h);
				*C = f2 /d;
				ਜੇ (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} ਨਹੀਂਤਾਂ {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			ਰੋਕੋ;
		}
	}
    } ਜਦੋਂ (0);

	FUNCTION_PROFILE_ਸਿਰਾ(4, 4, 4);
	IDEBUG_ਸਿਰਾ;
}
	
