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
શૂન્ય CNAME(શૂન્ય *VDA, શૂન્ય *VDB, FLOAT *C, શૂન્ય *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#endif /* CBLAS */

#ifdef DOUBLE
  long દ્વિગુણ safmin = DBL_MIN;
  long દ્વિગુણ rtmin = sqrt(DBL_MIN/DBL_EPSILON);
#નહીંતર
  long દ્વિગુણ safmin = FLT_MIN;
  long દ્વિગુણ rtmin = sqrt(FLT_MIN/FLT_EPSILON);
#endif


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  long દ્વિગુણ d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
  long દ્વિગુણ safmax = 1./safmin;
  safmax *= 2.;
  કરો {
	જો (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		તોડો;
	}

#જો defined DOUBLE
	long દ્વિગુણ rtmax = safmax /DBL_EPSILON;
#નહીંતર
	long દ્વિગુણ rtmax = safmax /FLT_EPSILON;
#endif

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	જો (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		જો (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			તોડો;
		} નહીંતર જો ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			તોડો;
		} નહીંતર {
			long દ્વિગુણ g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =sqrt(safmax/2.);
			જો (g1 > rtmin && g1 < rtmax) { // unscaled
				d = sqrt(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				તોડો;
			} નહીંતર { // scaled algorithm
				long દ્વિગુણ u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = sqrt ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				તોડો;
			}
		}

	} નહીંતર {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = sqrt(safmax / 4.);
		જો ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			long દ્વિગુણ h = ada + adb;
			દ્વિગુણ adahsq = sqrt(ada * h);
			જો (ada >= h *safmin) {
				*C = sqrt(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				જો ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} નહીંતર {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} નહીંતર {
				*C = ada / adahsq;
				જો (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} નહીંતર {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			તોડો;
		} નહીંતર {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			long દ્વિગુણ v,w,f2,g2,h;
			long દ્વિગુણ u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			જો (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} નહીંતર { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			જો ( f2 >= h * safmin) {
				*C = sqrt ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				જો ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=sqrt(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} નહીંતર {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} નહીંતર { // intermediates might overflow
				d = sqrt ( f2 * h);
				*C = f2 /d;
				જો (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} નહીંતર {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			તોડો;
		}
	}
    } જ્યારે (0);

	FUNCTION_PROFILE_END(4, 4, 4);
	IDEBUG_END;
}
	
