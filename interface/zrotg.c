#అన్తర్భావ <math.h>
#అన్తర్భావ <float.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు


OPENBLAS_EXPORT
#నిర్వచితం_కాకపోతే CBLAS
శూన్య NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#లేకపోతే
శూన్య CNAME(శూన్య *VDA, శూన్య *VDB, FLOAT *C, శూన్య *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#అయితే_ముగింపు /* CBLAS */

#నిర్వచితమైతే DOUBLE
  దీర్ఘ ద్విగుణ safmin = DBL_MIN;
  దీర్ఘ ద్విగుణ rtmin = వర్గమూలం(DBL_MIN/DBL_EPSILON);
#లేకపోతే
  దీర్ఘ ద్విగుణ safmin = FLT_MIN;
  దీర్ఘ ద్విగుణ rtmin = వర్గమూలం(FLT_MIN/FLT_EPSILON);
#అయితే_ముగింపు


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  దీర్ఘ ద్విగుణ d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_ప్రారంభం;

  FUNCTION_PROFILE_ప్రారంభం();
  దీర్ఘ ద్విగుణ safmax = 1./safmin;
  safmax *= 2.;
  చేయి {
	ఒకవేళ (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		ఆపు;
	}

#ఒకవేళ defined DOUBLE
	దీర్ఘ ద్విగుణ rtmax = safmax /DBL_EPSILON;
#లేకపోతే
	దీర్ఘ ద్విగుణ rtmax = safmax /FLT_EPSILON;
#అయితే_ముగింపు

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	ఒకవేళ (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		ఒకవేళ (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ఆపు;
		} లేకపోతే ఒకవేళ ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			ఆపు;
		} లేకపోతే {
			దీర్ఘ ద్విగుణ g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =వర్గమూలం(safmax/2.);
			ఒకవేళ (g1 > rtmin && g1 < rtmax) { // unscaled
				d = వర్గమూలం(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				ఆపు;
			} లేకపోతే { // scaled algorithm
				దీర్ఘ ద్విగుణ u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = వర్గమూలం ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				ఆపు;
			}
		}

	} లేకపోతే {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = వర్గమూలం(safmax / 4.);
		ఒకవేళ ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			దీర్ఘ ద్విగుణ h = ada + adb;
			ద్విగుణ adahsq = వర్గమూలం(ada * h);
			ఒకవేళ (ada >= h *safmin) {
				*C = వర్గమూలం(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				ఒకవేళ ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} లేకపోతే {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} లేకపోతే {
				*C = ada / adahsq;
				ఒకవేళ (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} లేకపోతే {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			ఆపు;
		} లేకపోతే {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			దీర్ఘ ద్విగుణ v,w,f2,g2,h;
			దీర్ఘ ద్విగుణ u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			ఒకవేళ (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} లేకపోతే { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			ఒకవేళ ( f2 >= h * safmin) {
				*C = వర్గమూలం ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				ఒకవేళ ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=వర్గమూలం(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} లేకపోతే {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} లేకపోతే { // intermediates might overflow
				d = వర్గమూలం ( f2 * h);
				*C = f2 /d;
				ఒకవేళ (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} లేకపోతే {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			ఆపు;
		}
	}
    } అయితే (0);

	FUNCTION_PROFILE_చివర(4, 4, 4);
	IDEBUG_చివర;
}
	
