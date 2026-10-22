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
#تعریف_کرو ERROR_NAME "XTBMV "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZTBMV "
#ورنہ
#تعریف_کرو ERROR_NAME "CTBMV "
#اگر_ختم

ساکن صحیح (*tbmv[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, خالی *) = {
#اگر_تعریف_شدہ XDOUBLE
  xtbmv_NUU, xtbmv_NUN, xtbmv_NLU, xtbmv_NLN,
  xtbmv_TUU, xtbmv_TUN, xtbmv_TLU, xtbmv_TLN,
  xtbmv_RUU, xtbmv_تسلسل, xtbmv_RLU, xtbmv_RLN,
  xtbmv_CUU, xtbmv_CUN, xtbmv_CLU, xtbmv_CLN,
#وگرنہ_اگر defined(DOUBLE)
  ztbmv_NUU, ztbmv_NUN, ztbmv_NLU, ztbmv_NLN,
  ztbmv_TUU, ztbmv_TUN, ztbmv_TLU, ztbmv_TLN,
  ztbmv_RUU, ztbmv_تسلسل, ztbmv_RLU, ztbmv_RLN,
  ztbmv_CUU, ztbmv_CUN, ztbmv_CLU, ztbmv_CLN,
#ورنہ
  ctbmv_NUU, ctbmv_NUN, ctbmv_NLU, ctbmv_NLN,
  ctbmv_TUU, ctbmv_TUN, ctbmv_TLU, ctbmv_TLN,
  ctbmv_RUU, ctbmv_تسلسل, ctbmv_RLU, ctbmv_RLN,
  ctbmv_CUU, ctbmv_CUN, ctbmv_CLU, ctbmv_CLN,
#اگر_ختم
};

#اگر_تعریف_شدہ SMP
ساکن صحیح (*tbmv_thread[])(BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT *, BLASLONG, FLOAT *, صحیح) = {
#اگر_تعریف_شدہ XDOUBLE
  xtbmv_thread_NUU, xtbmv_thread_NUN, xtbmv_thread_NLU, xtbmv_thread_NLN,
  xtbmv_thread_TUU, xtbmv_thread_TUN, xtbmv_thread_TLU, xtbmv_thread_TLN,
  xtbmv_thread_RUU, xtbmv_thread_تسلسل, xtbmv_thread_RLU, xtbmv_thread_RLN,
  xtbmv_thread_CUU, xtbmv_thread_CUN, xtbmv_thread_CLU, xtbmv_thread_CLN,
#وگرنہ_اگر defined(DOUBLE)
  ztbmv_thread_NUU, ztbmv_thread_NUN, ztbmv_thread_NLU, ztbmv_thread_NLN,
  ztbmv_thread_TUU, ztbmv_thread_TUN, ztbmv_thread_TLU, ztbmv_thread_TLN,
  ztbmv_thread_RUU, ztbmv_thread_تسلسل, ztbmv_thread_RLU, ztbmv_thread_RLN,
  ztbmv_thread_CUU, ztbmv_thread_CUN, ztbmv_thread_CLU, ztbmv_thread_CLN,
#ورنہ
  ctbmv_thread_NUU, ctbmv_thread_NUN, ctbmv_thread_NLU, ctbmv_thread_NLN,
  ctbmv_thread_TUU, ctbmv_thread_TUN, ctbmv_thread_TLU, ctbmv_thread_TLN,
  ctbmv_thread_RUU, ctbmv_thread_تسلسل, ctbmv_thread_RLU, ctbmv_thread_RLN,
  ctbmv_thread_CUU, ctbmv_thread_CUN, ctbmv_thread_CLU, ctbmv_thread_CLN,
#اگر_ختم
};
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *UPLO, حرف *TRANS, حرف *DIAG,
	 blasint *N, blasint *K,
	 FLOAT *a, blasint *LDA, FLOAT *x, blasint *INCX){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;
  حرف diag_arg  = *DIAG;

  blasint n    = *N;
  blasint k    = *K;
  blasint lda  = *LDA;
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
  اگر (trans_arg == 'R') trans = 2;
  اگر (trans_arg == 'C') trans = 3;

  اگر (diag_arg  == 'U') unit  = 0;
  اگر (diag_arg  == 'N') unit  = 1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;

  info = 0;

  اگر (incx == 0)          info =  9;
  اگر (lda < k + 1)        info =  7;
  اگر (k < 0)              info =  5;
  اگر (n < 0)              info =  4;
  اگر (unit  < 0)          info =  3;
  اگر (trans < 0)          info =  2;
  اگر (uplo  < 0)          info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_UPLO Uplo,
	   شمار CBLAS_TRANSPOSE TransA, شمار CBLAS_DIAG Diag,
	   blasint n, blasint k, خالی  *va, blasint lda, خالی  *vx, blasint incx) {

  FLOAT *a = (FLOAT*) va;
  FLOAT *x = (FLOAT*) vx;

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
    اگر (TransA == CblasConjNoTrans) trans = 2;
    اگر (TransA == CblasConjTrans)   trans = 3;

    اگر (Diag == CblasUnit)          unit  = 0;
    اگر (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    اگر (incx == 0)          info =  9;
    اگر (lda < k + 1)        info =  7;
    اگر (k < 0)              info =  5;
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
    اگر (TransA == CblasConjNoTrans) trans = 3;
    اگر (TransA == CblasConjTrans)   trans = 2;

    اگر (Diag == CblasUnit)          unit  = 0;
    اگر (Diag == CblasNonUnit)       unit  = 1;

    info = -1;

    اگر (incx == 0)          info =  9;
    اگر (lda < k + 1)        info =  7;
    اگر (k < 0)              info =  5;
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

  اگر (incx < 0 ) x -= (n - 1) * incx * 2;

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

  (tbmv[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, بفر);

#اگر_تعریف_شدہ SMP
  } ورنہ {

    (tbmv_thread[(trans<<2) | (uplo<<1) | unit])(n, k, a, lda, x, incx, بفر, nthreads);

  }
#اگر_ختم

  blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(4, n * k / 2 + n, n * k);

  IDEBUG_انتہا;

  واپس;
}
