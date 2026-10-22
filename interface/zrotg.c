#ଅନ୍ତର୍ଭାଵ <math.h>
#ଅନ୍ତର୍ଭାଵ <float.h>
#ଅନ୍ତର୍ଭାଵ "common.h"
#ifdef FUNCTION_PROFILE
#ଅନ୍ତର୍ଭାଵ "functable.h"
#endif


OPENBLAS_EXPORT
#ifndef CBLAS
ଶୂନ୍ୟ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ଅନ୍ୟଥା
ଶୂନ୍ୟ CNAME(ଶୂନ୍ୟ *VDA, ଶୂନ୍ୟ *VDB, FLOAT *C, ଶୂନ୍ୟ *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#endif /* CBLAS */

#ifdef DOUBLE
  long ଦ୍ୱିଗୁଣ safmin = DBL_MIN;
  long ଦ୍ୱିଗୁଣ rtmin = sqrt(DBL_MIN/DBL_EPSILON);
#ଅନ୍ୟଥା
  long ଦ୍ୱିଗୁଣ safmin = FLT_MIN;
  long ଦ୍ୱିଗୁଣ rtmin = sqrt(FLT_MIN/FLT_EPSILON);
#endif


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  long ଦ୍ୱିଗୁଣ d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
  long ଦ୍ୱିଗୁଣ safmax = 1./safmin;
  safmax *= 2.;
  କର {
	ଯଦି (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		ଭାଙ୍ଗ;
	}

#ଯଦି defined DOUBLE
	long ଦ୍ୱିଗୁଣ rtmax = safmax /DBL_EPSILON;
#ଅନ୍ୟଥା
	long ଦ୍ୱିଗୁଣ rtmax = safmax /FLT_EPSILON;
#endif

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	ଯଦି (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		ଯଦି (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ଭାଙ୍ଗ;
		} ଅନ୍ୟଥା ଯଦି ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ଭାଙ୍ଗ;
		} ଅନ୍ୟଥା {
			long ଦ୍ୱିଗୁଣ g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =sqrt(safmax/2.);
			ଯଦି (g1 > rtmin && g1 < rtmax) { // unscaled
				d = sqrt(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				ଭାଙ୍ଗ;
			} ଅନ୍ୟଥା { // scaled algorithm
				long ଦ୍ୱିଗୁଣ u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = sqrt ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				ଭାଙ୍ଗ;
			}
		}

	} ଅନ୍ୟଥା {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = sqrt(safmax / 4.);
		ଯଦି ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			long ଦ୍ୱିଗୁଣ h = ada + adb;
			ଦ୍ୱିଗୁଣ adahsq = sqrt(ada * h);
			ଯଦି (ada >= h *safmin) {
				*C = sqrt(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				ଯଦି ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} ଅନ୍ୟଥା {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} ଅନ୍ୟଥା {
				*C = ada / adahsq;
				ଯଦି (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} ଅନ୍ୟଥା {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			ଭାଙ୍ଗ;
		} ଅନ୍ୟଥା {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			long ଦ୍ୱିଗୁଣ v,w,f2,g2,h;
			long ଦ୍ୱିଗୁଣ u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			ଯଦି (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} ଅନ୍ୟଥା { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			ଯଦି ( f2 >= h * safmin) {
				*C = sqrt ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				ଯଦି ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=sqrt(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} ଅନ୍ୟଥା {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} ଅନ୍ୟଥା { // intermediates might overflow
				d = sqrt ( f2 * h);
				*C = f2 /d;
				ଯଦି (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} ଅନ୍ୟଥା {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			ଭାଙ୍ଗ;
		}
	}
    } ଯେତେବେଳେ (0);

	FUNCTION_PROFILE_END(4, 4, 4);
	IDEBUG_END;
}
	
