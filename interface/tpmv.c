/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#شامل <stdio.h>
#شامل <ctype.h>
#شامل "common.h"
#اگر_تعریف_شدہ FUNCTION_PROFILE
#شامل "functable.h"
#اگر_ختم

#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QTPMV "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DTPMV "
#ورنہ
#تعریف_کرو ERROR_NAME "STPMV "
#اگر_ختم

ساکن صحیح (*tpmv[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, خالی *) = {
#اگر_تعریف_شدہ XDOUBLE
  qtpmv_NUU, qtpmv_NUN, qtpmv_NLU, qtpmv_NLN,
  qtpmv_TUU, qtpmv_TUN, qtpmv_TLU, qtpmv_TLN,
#وگرنہ_اگر defined(DOUBLE)
  dtpmv_NUU, dtpmv_NUN, dtpmv_NLU, dtpmv_NLN,
  dtpmv_TUU, dtpmv_TUN, dtpmv_TLU, dtpmv_TLN,
#ورنہ
  stpmv_NUU, stpmv_NUN, stpmv_NLU, stpmv_NLN,
  stpmv_TUU, stpmv_TUN, stpmv_TLU, stpmv_TLN,
#اگر_ختم
};

#اگر_تعریف_شدہ SMP
ساکن صحیح (*tpmv_thread[])(BLASLONG, FLOAT *, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#اگر_تعریف_شدہ XDOUBLE
  qtpmv_thread_NUU, qtpmv_thread_NUN, qtpmv_thread_NLU, qtpmv_thread_NLN,
  qtpmv_thread_TUU, qtpmv_thread_TUN, qtpmv_thread_TLU, qtpmv_thread_TLN,
#وگرنہ_اگر defined(DOUBLE)
  dtpmv_thread_NUU, dtpmv_thread_NUN, dtpmv_thread_NLU, dtpmv_thread_NLN,
  dtpmv_thread_TUU, dtpmv_thread_TUN, dtpmv_thread_TLU, dtpmv_thread_TLN,
#ورنہ
  stpmv_thread_NUU, stpmv_thread_NUN, stpmv_thread_NLU, stpmv_thread_NLN,
  stpmv_thread_TUU, stpmv_thread_TUN, stpmv_thread_TLU, stpmv_thread_TLN,
#اگر_ختم
};
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *UPLO, حرف *TRANS, حرف *DIAG,
	   blasint *N, FLOAT *a, FLOAT *x, blasint *INCX){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blasint n    = *N;
  blasint incx = *INCX;

  blasint info;
  صحیح uplo;
  صحیح unit;
  صحیح trans;
  FLOAT *بفر;
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

  PRINT_DEBUG_NAME;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);
  TOUPPER(diag_arg);

  trans = -1;
  unit  = -1;
  uplo  = -1;

  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
  اگر (trans_arg == 'R') trans = 0;
  اگر (trans_arg == 'C') trans = 1;

  اگر (diag_arg  == 'U') unit  = 0;
  اگر (diag_arg  == 'N') unit  = 1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  اگر (incx == 0)          info =  7;
  اگر (n < 0)              info =  4;
  اگر (unit  < 0)          info =  3;
  اگر (trans < 0)          info =  2;
  اگر (uplo  < 0)          info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ
#اگر_غیر_تعریف_شدہ COMPLEX
خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_UPLO Uplo,
	   شمار CBLAS_TRANSPOSE TransA, شمار CBLAS_DIAG Diag,
	   blasint n, FLOAT  *a, FLOAT  *x, blasint incx) {
#ورنہ
خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_UPLO Uplo,
	   شمار CBLAS_TRANSPOSE TransA, شمار CBLAS_DIAG Diag,
	   blasint n, خالی  *va, خالی  *vx, blasint incx) {
  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;
#اگر_ختم

  صحیح trans, uplo, unit;
  blasint info;
  FLOAT *بفر;
#اگر_تعریف_شدہ SMP
  صحیح nthreads;
#اگر_ختم

  PRINT_DEBUG_CNAME;

  unit  = -1;
  uplo  = -1;
  trans = -1;
  info  =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
    اگر (Uplo == Cblas_بالائی)         uplo  = 0;
    اگر (Uplo == Cblas_زیریں)         uplo  = 1;

    اگر (TransA == CblasNoTrans)     trans = 0;
    اگر (TransA == CblasTrans)       trans = 1;
    اگر (TransA == CblasConjNoTrans) trans = 0;
    اگر (TransA == CblasConjTrans)   trans = 1;

    اگر (Diag == CblasUnit)          unit  = 0;
    اگر (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    اگر (incx == 0)          info =  7;
    اگر (n < 0)              info =  4;
    اگر (unit  < 0)          info =  3;
    اگر (trans < 0)          info =  2;
    اگر (uplo  < 0)          info =  1;
  }

  اگر (ترتیب == Cblas_صف_Major) {
    اگر (Uplo == Cblas_بالائی)         uplo  = 1;
    اگر (Uplo == Cblas_زیریں)         uplo  = 0;

    اگر (TransA == CblasNoTrans)     trans = 1;
    اگر (TransA == CblasTrans)       trans = 0;
    اگر (TransA == CblasConjNoTrans) trans = 1;
    اگر (TransA == CblasConjTrans)   trans = 0;

    اگر (Diag == CblasUnit)          unit  = 0;
    اگر (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    اگر (incx == 0)          info =  7;
    اگر (n < 0)              info =  4;
    اگر (unit  < 0)          info =  3;
    اگر (trans < 0)          info =  2;
    اگر (uplo  < 0)          info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_ختم

  اگر (n == 0) واپس;

  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  اگر (incx < 0 ) x -= (n - 1) * incx;

  بفر = (FLOAT *)blas_حافظہ_alloc(1);
  اگر (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#اگر_تعریف_شدہ SMP
  nthreads = عدد_cpu_avail(2);

  اگر (n < 50) nthreads = 1;
  اگر (nthreads > 2 && n < 500) nthreads = 2;

  اگر (nthreads == 1) {
#اگر_ختم

  (tpmv[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, بفر);

#اگر_تعریف_شدہ SMP
  } ورنہ {

    (tpmv_thread[(trans<<2) | (uplo<<1) | unit])(n, a, x, incx, بفر, nthreads);

  }
#اگر_ختم

  blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(1, n * n / 2 + n, n * n);

  IDEBUG_انتہا;

  واپس;
}
