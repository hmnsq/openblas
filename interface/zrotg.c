#include <math.h>
#include <float.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif


OPENBLAS_EXPORT
#ifndef CBLAS
void NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#ບໍ່ດັ່ງນັ້ນ
void CNAME(void *VDA, void *VDB, FLOAT *C, void *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#endif /* CBLAS */

#ifdef DOUBLE
  long dviguu_na safmin = DBL_MIN;
  long dviguu_na rtmin = sqrt(DBL_MIN/DBL_EPSILON);
#ບໍ່ດັ່ງນັ້ນ
  long dviguu_na safmin = FLT_MIN;
  long dviguu_na rtmin = sqrt(FLT_MIN/FLT_EPSILON);
#endif


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  long dviguu_na d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
  long dviguu_na safmax = 1./safmin;
  safmax *= 2.;
  do {
	ຖ້າ (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		bha_nga;
	}

#ຖ້າ defined DOUBLE
	long dviguu_na rtmax = safmax /DBL_EPSILON;
#ບໍ່ດັ່ງນັ້ນ
	long dviguu_na rtmax = safmax /FLT_EPSILON;
#endif

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	ຖ້າ (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		ຖ້າ (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			bha_nga;
		} ບໍ່ດັ່ງນັ້ນ ຖ້າ ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			bha_nga;
		} ບໍ່ດັ່ງນັ້ນ {
			long dviguu_na g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =sqrt(safmax/2.);
			ຖ້າ (g1 > rtmin && g1 < rtmax) { // unscaled
				d = sqrt(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				bha_nga;
			} ບໍ່ດັ່ງນັ້ນ { // scaled algorithm
				long dviguu_na u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = sqrt ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				bha_nga;
			}
		}

	} ບໍ່ດັ່ງນັ້ນ {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = sqrt(safmax / 4.);
		ຖ້າ ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			long dviguu_na h = ada + adb;
			dviguu_na adahsq = sqrt(ada * h);
			ຖ້າ (ada >= h *safmin) {
				*C = sqrt(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				ຖ້າ ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} ບໍ່ດັ່ງນັ້ນ {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} ບໍ່ດັ່ງນັ້ນ {
				*C = ada / adahsq;
				ຖ້າ (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} ບໍ່ດັ່ງນັ້ນ {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			bha_nga;
		} ບໍ່ດັ່ງນັ້ນ {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			long dviguu_na v,w,f2,g2,h;
			long dviguu_na u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			ຖ້າ (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} ບໍ່ດັ່ງນັ້ນ { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			ຖ້າ ( f2 >= h * safmin) {
				*C = sqrt ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				ຖ້າ ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=sqrt(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} ບໍ່ດັ່ງນັ້ນ {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} ບໍ່ດັ່ງນັ້ນ { // intermediates might overflow
				d = sqrt ( f2 * h);
				*C = f2 /d;
				ຖ້າ (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} ບໍ່ດັ່ງນັ້ນ {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			bha_nga;
		}
	}
    } ໃນຂະນະທີ່ (0);

	FUNCTION_PROFILE_END(4, 4, 4);
	IDEBUG_END;
}
	
