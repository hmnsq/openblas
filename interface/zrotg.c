#ಅನ್ತರ್ಭಾವ <math.h>
#ಅನ್ತರ್ಭಾವ <float.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#ಆದರೆ_ಅಂತ್ಯ


OPENBLAS_EXPORT
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
ಶೂನ್ಯ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ಇಲ್ಲದಿದ್ದರೆ
ಶೂನ್ಯ CNAME(ಶೂನ್ಯ *VDA, ಶೂನ್ಯ *VDB, FLOAT *C, ಶೂನ್ಯ *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#ಆದರೆ_ಅಂತ್ಯ /* CBLAS */

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ DOUBLE
  ದೀರ್ಘ ದ್ವಿಗುಣ safmin = DBL_MIN;
  ದೀರ್ಘ ದ್ವಿಗುಣ rtmin = ವರ್ಗಮೂಲ(DBL_MIN/DBL_EPSILON);
#ಇಲ್ಲದಿದ್ದರೆ
  ದೀರ್ಘ ದ್ವಿಗುಣ safmin = FLT_MIN;
  ದೀರ್ಘ ದ್ವಿಗುಣ rtmin = ವರ್ಗಮೂಲ(FLT_MIN/FLT_EPSILON);
#ಆದರೆ_ಅಂತ್ಯ


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  ದೀರ್ಘ ದ್ವಿಗುಣ d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_ಆರಂಭ;

  FUNCTION_PROFILE_ಆರಂಭ();
  ದೀರ್ಘ ದ್ವಿಗುಣ safmax = 1./safmin;
  safmax *= 2.;
  ಮಾಡು {
	ಆದರೆ (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		ನಿಲ್ಲಿಸು;
	}

#ಆದರೆ defined DOUBLE
	ದೀರ್ಘ ದ್ವಿಗುಣ rtmax = safmax /DBL_EPSILON;
#ಇಲ್ಲದಿದ್ದರೆ
	ದೀರ್ಘ ದ್ವಿಗುಣ rtmax = safmax /FLT_EPSILON;
#ಆದರೆ_ಅಂತ್ಯ

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	ಆದರೆ (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		ಆದರೆ (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ನಿಲ್ಲಿಸು;
		} ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ನಿಲ್ಲಿಸು;
		} ಇಲ್ಲದಿದ್ದರೆ {
			ದೀರ್ಘ ದ್ವಿಗುಣ g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =ವರ್ಗಮೂಲ(safmax/2.);
			ಆದರೆ (g1 > rtmin && g1 < rtmax) { // unscaled
				d = ವರ್ಗಮೂಲ(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				ನಿಲ್ಲಿಸು;
			} ಇಲ್ಲದಿದ್ದರೆ { // scaled algorithm
				ದೀರ್ಘ ದ್ವಿಗುಣ u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = ವರ್ಗಮೂಲ ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				ನಿಲ್ಲಿಸು;
			}
		}

	} ಇಲ್ಲದಿದ್ದರೆ {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = ವರ್ಗಮೂಲ(safmax / 4.);
		ಆದರೆ ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			ದೀರ್ಘ ದ್ವಿಗುಣ h = ada + adb;
			ದ್ವಿಗುಣ adahsq = ವರ್ಗಮೂಲ(ada * h);
			ಆದರೆ (ada >= h *safmin) {
				*C = ವರ್ಗಮೂಲ(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				ಆದರೆ ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} ಇಲ್ಲದಿದ್ದರೆ {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} ಇಲ್ಲದಿದ್ದರೆ {
				*C = ada / adahsq;
				ಆದರೆ (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} ಇಲ್ಲದಿದ್ದರೆ {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			ನಿಲ್ಲಿಸು;
		} ಇಲ್ಲದಿದ್ದರೆ {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			ದೀರ್ಘ ದ್ವಿಗುಣ v,w,f2,g2,h;
			ದೀರ್ಘ ದ್ವಿಗುಣ u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			ಆದರೆ (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} ಇಲ್ಲದಿದ್ದರೆ { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			ಆದರೆ ( f2 >= h * safmin) {
				*C = ವರ್ಗಮೂಲ ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				ಆದರೆ ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=ವರ್ಗಮೂಲ(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} ಇಲ್ಲದಿದ್ದರೆ {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} ಇಲ್ಲದಿದ್ದರೆ { // intermediates might overflow
				d = ವರ್ಗಮೂಲ ( f2 * h);
				*C = f2 /d;
				ಆದರೆ (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} ಇಲ್ಲದಿದ್ದರೆ {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			ನಿಲ್ಲಿಸು;
		}
	}
    } ಆಗ (0);

	FUNCTION_PROFILE_ಕೊನೆ(4, 4, 4);
	IDEBUG_ಕೊನೆ;
}
	
