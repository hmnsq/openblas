#ضمن <math.h>
#ضمن <float.h>
#ضمن "common.h"
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت


OPENBLAS_EXPORT
#جے_بنا_تعریف CBLAS
فراغ NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#وإلا
فراغ CNAME(فراغ *VDA, فراغ *VDB, FLOAT *C, فراغ *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#جے_انت /* CBLAS */

#جے_تعریف_شدہ DOUBLE
  لما مزدوج safmin = DBL_MIN;
  لما مزدوج rtmin = جذر(DBL_MIN/DBL_EPSILON);
#وإلا
  لما مزدوج safmin = FLT_MIN;
  لما مزدوج rtmin = جذر(FLT_MIN/FLT_EPSILON);
#جے_انت


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  لما مزدوج d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();
  لما مزدوج safmax = 1./safmin;
  safmax *= 2.;
  افعل {
	إذا (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		اكسر;
	}

#إذا defined DOUBLE
	لما مزدوج rtmax = safmax /DBL_EPSILON;
#وإلا
	لما مزدوج rtmax = safmax /FLT_EPSILON;
#جے_انت

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	إذا (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		إذا (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			اكسر;
		} وإلا إذا ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			اكسر;
		} وإلا {
			لما مزدوج g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =جذر(safmax/2.);
			إذا (g1 > rtmin && g1 < rtmax) { // unscaled
				d = جذر(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				اكسر;
			} وإلا { // scaled algorithm
				لما مزدوج u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = جذر ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				اكسر;
			}
		}

	} وإلا {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = جذر(safmax / 4.);
		إذا ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			لما مزدوج h = ada + adb;
			مزدوج adahsq = جذر(ada * h);
			إذا (ada >= h *safmin) {
				*C = جذر(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				إذا ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} وإلا {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} وإلا {
				*C = ada / adahsq;
				إذا (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} وإلا {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			اكسر;
		} وإلا {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			لما مزدوج v,w,f2,g2,h;
			لما مزدوج u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			إذا (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} وإلا { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			إذا ( f2 >= h * safmin) {
				*C = جذر ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				إذا ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=جذر(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} وإلا {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} وإلا { // intermediates might overflow
				d = جذر ( f2 * h);
				*C = f2 /d;
				إذا (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} وإلا {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			اكسر;
		}
	}
    } طالما (0);

	FUNCTION_PROFILE_سرا(4, 4, 4);
	IDEBUG_سرا;
}
	
