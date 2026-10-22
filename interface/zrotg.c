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
خالی CNAME(خالی *VDA, خالی *VDB, FLOAT *C, خالی *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#اگر_ختم /* CBLAS */

#اگر_تعریف_شدہ DOUBLE
  طویل دوہرا safmin = DBL_MIN;
  طویل دوہرا rtmin = جذر(DBL_MIN/DBL_EPSILON);
#ورنہ
  طویل دوہرا safmin = FLT_MIN;
  طویل دوہرا rtmin = جذر(FLT_MIN/FLT_EPSILON);
#اگر_ختم


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  طویل دوہرا d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();
  طویل دوہرا safmax = 1./safmin;
  safmax *= 2.;
  کرو {
	اگر (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		توڑو;
	}

#اگر defined DOUBLE
	طویل دوہرا rtmax = safmax /DBL_EPSILON;
#ورنہ
	طویل دوہرا rtmax = safmax /FLT_EPSILON;
#اگر_ختم

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	اگر (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		اگر (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			توڑو;
		} ورنہ اگر ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			توڑو;
		} ورنہ {
			طویل دوہرا g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =جذر(safmax/2.);
			اگر (g1 > rtmin && g1 < rtmax) { // unscaled
				d = جذر(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				توڑو;
			} ورنہ { // scaled algorithm
				طویل دوہرا u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = جذر ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				توڑو;
			}
		}

	} ورنہ {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = جذر(safmax / 4.);
		اگر ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			طویل دوہرا h = ada + adb;
			دوہرا adahsq = جذر(ada * h);
			اگر (ada >= h *safmin) {
				*C = جذر(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				اگر ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} ورنہ {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} ورنہ {
				*C = ada / adahsq;
				اگر (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} ورنہ {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			توڑو;
		} ورنہ {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			طویل دوہرا v,w,f2,g2,h;
			طویل دوہرا u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			اگر (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} ورنہ { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			اگر ( f2 >= h * safmin) {
				*C = جذر ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				اگر ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=جذر(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} ورنہ {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} ورنہ { // intermediates might overflow
				d = جذر ( f2 * h);
				*C = f2 /d;
				اگر (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} ورنہ {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			توڑو;
		}
	}
    } جب (0);

	FUNCTION_PROFILE_انتہا(4, 4, 4);
	IDEBUG_انتہا;
}
	
