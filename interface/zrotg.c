#अन्तर्भाव <math.h>
#अन्तर्भाव <float.h>
#अन्तर्भाव "common.h"
#ifdef FUNCTION_PROFILE
#अन्तर्भाव "functable.h"
#endif


OPENBLAS_EXPORT
#ifndef CBLAS
शून्य NAME(FLOAT *DA, FLOAT *DB, FLOAT *C, FLOAT *S){

#अन्यथा
शून्य CNAME(शून्य *VDA, शून्य *VDB, FLOAT *C, शून्य *VS) {
    FLOAT *DA = (FLOAT*) VDA;
    FLOAT *DB = (FLOAT*) VDB;
    FLOAT *S  = (FLOAT*) VS;
#endif /* CBLAS */

#ifdef DOUBLE
  दीर्घ द्विगुण safmin = द्विग_न्यून;
  दीर्घ द्विगुण rtmin = वर्ग_मूल(द्विग_न्यून/द्विग_अंतर);
#अन्यथा
  दीर्घ द्विगुण safmin = भग्न_न्यून;
  दीर्घ द्विगुण rtmin = वर्ग_मूल(भग्न_न्यून/भग्न_अंतर);
#endif


  FLOAT da_r = *(DA+0);
  FLOAT da_i = *(DA+1);
  FLOAT db_r = *(DB+0);
  FLOAT db_i = *(DB+1);
  //long double r;
  FLOAT S1[2];
  FLOAT R[2];
  दीर्घ द्विगुण d;

  FLOAT ada =  da_r * da_r + da_i * da_i; 
  FLOAT adb =  db_r * db_r + db_i * db_i; 

  PRINT_DEBUG_NAME;

  IDEBUG_START;

  FUNCTION_PROFILE_START();
  दीर्घ द्विगुण safmax = 1./safmin;
  safmax *= 2.;
  करा {
	जर (db_r == ZERO && db_i == ZERO) {
		*C        = ONE;
		*(S  + 0) = ZERO;
		*(S  + 1) = ZERO;
		थांबा;
	}

#जर defined DOUBLE
	दीर्घ द्विगुण rtmax = safmax /द्विग_अंतर;
#अन्यथा
	दीर्घ द्विगुण rtmax = safmax /भग्न_अंतर;
#endif

	*(S1 + 0) = *(DB + 0);
	*(S1 + 1) = *(DB + 1) *-1;

	जर (da_r == ZERO && da_i == ZERO) {
		*C = ZERO;
		जर (db_r == ZERO) {
			(*DA) = fabsl(db_i);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			थांबा;
		} अन्यथा जर ( db_i == ZERO) {
			*DA = fabsl(db_r);
			*S = *S1 /(*DA);
			*(S+1) = *(S1+1) /(*DA);
			थांबा;
		} अन्यथा {
			दीर्घ द्विगुण g1 = MAX( fabsl(db_r), fabsl(db_i));
			rtmax =वर्ग_मूल(safmax/2.);
			जर (g1 > rtmin && g1 < rtmax) { // unscaled
				d = वर्ग_मूल(adb);
				*S = *S1 /d;
				*(S+1) = *(S1+1) /d;
				*DA = d ;
				*(DA+1) = ZERO;
				थांबा;
			} अन्यथा { // scaled algorithm
				दीर्घ द्विगुण u = MIN ( safmax, MAX ( safmin, g1));
				FLOAT gs_r = db_r/u;
				FLOAT gs_i = db_i/u;
				d = वर्ग_मूल ( gs_r*gs_r + gs_i*gs_i);
				*S = gs_r / d;
				*(S + 1) = (gs_i * -1) / d;
				*DA = d * u;
				*(DA+1) = ZERO;
				थांबा;
			}
		}

	} अन्यथा {

		FLOAT f1 = MAX ( fabsl(da_r), fabsl(da_i));
		FLOAT g1 = MAX ( fabsl(db_r), fabsl(db_i));
		rtmax = वर्ग_मूल(safmax / 4.);
		जर ( f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) { //unscaled
			दीर्घ द्विगुण h = ada + adb;
			द्विगुण adahsq = वर्ग_मूल(ada * h);
			जर (ada >= h *safmin) {
				*C = वर्ग_मूल(ada/h);
				*R = *DA / *C;
				*(R+1) = *(DA+1) / *C;
				rtmax *= 2.;
				जर ( ada > rtmin && h < rtmax) { // no risk of intermediate overflow
					*S = *S1 * (*DA / adahsq) - *(S1+1)* (*(DA+1)/adahsq);
					*(S+1) = *S1 * (*(DA+1) / adahsq) + *(S1+1) * (*DA/adahsq);
				} अन्यथा {
					*S = *S1 * (*R/h) - *(S1+1) * (*(R+1)/h);
					*(S+1) = *S1 * (*(R+1)/h) + *(S1+1) * (*(R)/h);
				}
			} अन्यथा {
				*C = ada / adahsq;
				जर (*C >= safmin) {
					*R = *DA / *C;
					*(R+1) = *(DA+1) / *C;
				} अन्यथा {
					*R = *DA * (h / adahsq);
					*(R+1) = *(DA+1) * (h / adahsq);
				}
				*S = *S1 * (*DA) / adahsq;
				*(S+1) = *(S1+1) * (*DA) / adahsq;
			}
			*DA=*R;
			*(DA+1)=*(R+1);
			थांबा;
		} अन्यथा {	// scaled
			FLOAT fs_r, fs_i, gs_r, gs_i;
			दीर्घ द्विगुण v,w,f2,g2,h;
			दीर्घ द्विगुण u = MIN ( safmax, MAX ( safmin, MAX(f1,g1)));
			gs_r = db_r/u;
			gs_i = db_i/u;
			g2 = gs_r*gs_r + gs_i*gs_i;
			जर (f1 /u < rtmin) {
				v = MIN (safmax, MAX (safmin, f1));
				w = v / u;
				fs_r = *DA/ v;
				fs_i = *(DA+1) / v;
				f2 = fs_r*fs_r + fs_i*fs_i;
				h = f2 * w * w + g2;		
			} अन्यथा { // use same scaling for both
				w = 1.;
				fs_r = *DA/ u;
				fs_i = *(DA+1) / u;
				f2 =  fs_r*fs_r + fs_i*fs_i;
				h = f2 + g2;
			}
			जर ( f2 >= h * safmin) {
				*C = वर्ग_मूल ( f2 / h );
				*DA = fs_r / *C;
				*(DA+1) = fs_i / *C;
				rtmax *= 2;
				जर ( f2 > rtmin && h < rtmax) {
					FLOAT f2h=वर्ग_मूल(f2*h);
					*S = gs_r * (fs_r /f2h) + gs_i * (fs_i / f2h);
					*(S+1) = gs_r * (fs_i /f2h) - gs_i * -1. * (fs_r / f2h);
				} अन्यथा {
					*S = gs_r * (*DA/h) - gs_i * (*(DA+1) / h);
					*(S+1) = gs_r * (*(DA+1) /h) + gs_i * -1. * (*DA / h);
				}
			} अन्यथा { // intermediates might overflow
				d = वर्ग_मूल ( f2 * h);
				*C = f2 /d;
				जर (*C >= safmin) {
					*DA = fs_r / *C;
					*(DA+1) = fs_i / *C;
				} अन्यथा {
					*DA = fs_r * (h / d);
					*(DA+1) = fs_i / (h / d);
				}
				*S = gs_r * (fs_r /d) - gs_i * (fs_i / d);
				*(S+1) = gs_r * (fs_i /d) + gs_i * -1. * (fs_r / d);
			}
			*C *= w;
			*DA *= u;
			*(DA+1) *= u;
			थांबा;
		}
	}
    } पर्यंत (0);

	FUNCTION_PROFILE_END(4, 4, 4);
	IDEBUG_END;
}
	
