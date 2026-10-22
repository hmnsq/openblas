#അന്തര്ഭാവ <math.h>
#അന്തര്ഭാവ <float.h>
#അന്തര്ഭാവ "common.h"
#ifdef FUNCTION_PROFILE
#അന്തര്ഭാവ "functable.h"
#endif


OPENBLAS_EXPORT
#ifndef CBLAS
ശൂന്യം NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#അല്ലെങ്കിൽ
ശൂന്യം CNAME(ശൂന്യം *VDA, ശൂന്യം *VDB, FLOAT *C, ശൂന്യം *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#endif /* CBLAS */

#ifdef DOUBLE
  long ഇരട്ട safmin = DBL_MIN;
  long ഇരട്ട rtmin = sqrt(DBL_MIN/DBL_EPSILON);
#അല്ലെങ്കിൽ
  long ഇരട്ട safmin = FLT_MIN;
  long ഇരട്ട rtmin = sqrt(FLT_MIN/FLT_EPSILON);
#endif


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  long ഇരട്ട d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
  long ഇരട്ട safmax = 1./safmin;
  safmax *= 2.;
  ചെയ്യുക {
	എങ്കിൽ (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		നിർത്തുക;
	}

#എങ്കിൽ defined DOUBLE
	long ഇരട്ട rtmax = safmax /DBL_EPSILON;
#അല്ലെങ്കിൽ
	long ഇരട്ട rtmax = safmax /FLT_EPSILON;
#endif

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	എങ്കിൽ (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		എങ്കിൽ (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			നിർത്തുക;
		} അല്ലെങ്കിൽ എങ്കിൽ ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			നിർത്തുക;
		} അല്ലെങ്കിൽ {
			long ഇരട്ട g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =sqrt(safmax/2.);
			എങ്കിൽ (g1 > rtmin && g1 < rtmax) { // unscaled
				d = sqrt(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				നിർത്തുക;
			} അല്ലെങ്കിൽ { // scaled algorithm
				long ഇരട്ട u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = sqrt ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				നിർത്തുക;
			}
		}

	} അല്ലെങ്കിൽ {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = sqrt(safmax / 4.);
		എങ്കിൽ ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			long ഇരട്ട h = ada + adb;
			ഇരട്ട adahsq = sqrt(ada * h);
			എങ്കിൽ (ada >= h *safmin) {
				*C = sqrt(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				എങ്കിൽ ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} അല്ലെങ്കിൽ {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} അല്ലെങ്കിൽ {
				*C = ada / adahsq;
				എങ്കിൽ (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} അല്ലെങ്കിൽ {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			നിർത്തുക;
		} അല്ലെങ്കിൽ {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			long ഇരട്ട v,w,f2,g2,h;
			long ഇരട്ട u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			എങ്കിൽ (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} അല്ലെങ്കിൽ { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			എങ്കിൽ ( f2 >= h * safmin) {
				*C = sqrt ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				എങ്കിൽ ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=sqrt(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} അല്ലെങ്കിൽ {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} അല്ലെങ്കിൽ { // intermediates might overflow
				d = sqrt ( f2 * h);
				*C = f2 /d;
				എങ്കിൽ (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} അല്ലെങ്കിൽ {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			നിർത്തുക;
		}
	}
    } ആകുമ്പോൾ (0);

	FUNCTION_PROFILE_END(4, 4, 4);
	IDEBUG_END;
}
	
