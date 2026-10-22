#ಅನ್ತರ್ಭಾವ <math.h>
#ಅನ್ತರ್ಭಾವ <float.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif


OPENBLAS_EXPORT
#ifndef CBLAS
ಶೂನ್ಯ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ಅನ್ಯಥಾ
ಶೂನ್ಯ CNAME(ಶೂನ್ಯ *VDA, ಶೂನ್ಯ *VDB, FLOAT *C, ಶೂನ್ಯ *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#endif /* CBLAS */

#ifdef DOUBLE
  long ದ್ವಿಗುಣ safmin = DBL_MIN;
  long ದ್ವಿಗುಣ rtmin = sqrt(DBL_MIN/DBL_EPSILON);
#ಅನ್ಯಥಾ
  long ದ್ವಿಗುಣ safmin = FLT_MIN;
  long ದ್ವಿಗುಣ rtmin = sqrt(FLT_MIN/FLT_EPSILON);
#endif


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  long ದ್ವಿಗುಣ d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
  long ದ್ವಿಗುಣ safmax = 1./safmin;
  safmax *= 2.;
  ಕುರು {
	ಯದಿ (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		ಭಙ್ಗ;
	}

#ಯದಿ defined DOUBLE
	long ದ್ವಿಗುಣ rtmax = safmax /DBL_EPSILON;
#ಅನ್ಯಥಾ
	long ದ್ವಿಗುಣ rtmax = safmax /FLT_EPSILON;
#endif

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	ಯದಿ (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		ಯದಿ (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ಭಙ್ಗ;
		} ಅನ್ಯಥಾ ಯದಿ ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ಭಙ್ಗ;
		} ಅನ್ಯಥಾ {
			long ದ್ವಿಗುಣ g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =sqrt(safmax/2.);
			ಯದಿ (g1 > rtmin && g1 < rtmax) { // unscaled
				d = sqrt(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				ಭಙ್ಗ;
			} ಅನ್ಯಥಾ { // scaled algorithm
				long ದ್ವಿಗುಣ u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = sqrt ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				ಭಙ್ಗ;
			}
		}

	} ಅನ್ಯಥಾ {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = sqrt(safmax / 4.);
		ಯದಿ ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			long ದ್ವಿಗುಣ h = ada + adb;
			ದ್ವಿಗುಣ adahsq = sqrt(ada * h);
			ಯದಿ (ada >= h *safmin) {
				*C = sqrt(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				ಯದಿ ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} ಅನ್ಯಥಾ {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} ಅನ್ಯಥಾ {
				*C = ada / adahsq;
				ಯದಿ (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} ಅನ್ಯಥಾ {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			ಭಙ್ಗ;
		} ಅನ್ಯಥಾ {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			long ದ್ವಿಗುಣ v,w,f2,g2,h;
			long ದ್ವಿಗುಣ u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			ಯದಿ (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} ಅನ್ಯಥಾ { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			ಯದಿ ( f2 >= h * safmin) {
				*C = sqrt ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				ಯದಿ ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=sqrt(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} ಅನ್ಯಥಾ {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} ಅನ್ಯಥಾ { // intermediates might overflow
				d = sqrt ( f2 * h);
				*C = f2 /d;
				ಯದಿ (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} ಅನ್ಯಥಾ {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			ಭಙ್ಗ;
		}
	}
    } ಯಾವತ್ (0);

	FUNCTION_PROFILE_END(4, 4, 4);
	IDEBUG_END;
}
	
