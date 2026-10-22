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

#اگر_غیر_تعریف_شدہ COMPLEX
#تعریف_کرو SMP_حد_MIN 109944.
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QSYRK "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DSYRK "
#ورنہ
#تعریف_کرو ERROR_NAME "SSYRK "
#اگر_ختم
#ورنہ
#تعریف_کرو SMP_حد_MIN 14824.
#اگر_غیر_تعریف_شدہ HEMM
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XSYRK "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZSYRK "
#ورنہ
#تعریف_کرو ERROR_NAME "CSYRK "
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XHERK "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZHERK "
#ورنہ
#تعریف_کرو ERROR_NAME "CHERK "
#اگر_ختم
#اگر_ختم
#اگر_ختم

#اگر_غیر_تعریف_شدہ GEMM_MULTITHREAD_حد
#تعریف_کرو GEMM_MULTITHREAD_حد 4
#اگر_ختم

ساکن صحیح (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#اگر_غیر_تعریف_شدہ HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#اگر_ختم
#ورنہ
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#اگر defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#اگر_ختم
#اگر_ختم
};

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS

خالی NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t دلائل;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#اگر_تعریف_شدہ SMP
  دوہرا NNK;
#اگر_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم
#اگر_ختم

  blasint info;
  صحیح uplo;
  صحیح trans;
  صحیح nrowa;

  PRINT_DEBUG_NAME;

  دلائل.n = *N;
  دلائل.k = *K;

  دلائل.a = (خالی *)a;
  دلائل.c = (خالی *)c;

  دلائل.lda = *ldA;
  دلائل.ldc = *ldC;

  دلائل.alpha = (خالی *)alpha;
  دلائل.beta  = (خالی *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  اگر (uplo_arg  == 'U') uplo  = 0;
  اگر (uplo_arg  == 'L') uplo  = 1;


#اگر_غیر_تعریف_شدہ COMPLEX
  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
  اگر (trans_arg == 'C') trans = 1;
#ورنہ
#اگر_تعریف_شدہ HEMM
  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'C') trans = 1;
#ورنہ
  اگر (trans_arg == 'N') trans = 0;
  اگر (trans_arg == 'T') trans = 1;
#اگر_ختم

#اگر_ختم

  nrowa = دلائل.n;
  اگر (trans & 1) nrowa = دلائل.k;

  info = 0;

  اگر (دلائل.ldc < MAX(1,دلائل.n)) info = 10;
  اگر (دلائل.lda < MAX(1,nrowa))  info =  7;
  اگر (دلائل.k < 0)               info =  4;
  اگر (دلائل.n < 0)               info =  3;
  اگر (trans < 0)                info =  2;
  اگر (uplo  < 0)                info =  1;

  اگر (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_UPLO Uplo, شمار CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#اگر !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#ورنہ
	   خالی *valpha,
#اگر_ختم
#اگر !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#ورنہ
	   خالی *va, blasint lda,
#اگر_ختم
#اگر !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ورنہ
	   خالی *vbeta,
#اگر_ختم
#اگر !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#ورنہ
	   خالی *vc, blasint ldc) {
#اگر_ختم

#اگر_تعریف_شدہ COMPLEX
#اگر !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#اگر_ختم
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#اگر_ختم

  blas_arg_t دلائل;
  صحیح uplo, trans;
  blasint info, nrowa;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#اگر_تعریف_شدہ SMP
دوہرا NNK;

#اگر_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
  صحیح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
  صحیح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
  صحیح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم
#اگر_ختم

  PRINT_DEBUG_CNAME;

  دلائل.n = n;
  دلائل.k = k;

  دلائل.a = (خالی *)a;
  دلائل.c = (خالی *)c;

  دلائل.lda = lda;
  دلائل.ldc = ldc;

#اگر !defined(COMPLEX) || defined(HEMM)
  دلائل.alpha = (خالی *)&alpha;
  دلائل.beta  = (خالی *)&beta;
#ورنہ
  دلائل.alpha = (خالی *)alpha;
  دلائل.beta  = (خالی *)beta;
#اگر_ختم

  trans = -1;
  uplo  = -1;
  info  =  0;

  اگر (ترتیب == Cblas_ستون_Major) {
    اگر (Uplo == Cblas_بالائی) uplo  = 0;
    اگر (Uplo == Cblas_زیریں) uplo  = 1;

    اگر (Trans == CblasNoTrans)     trans = 0;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (Trans == CblasTrans)       trans = 1;
    اگر (Trans == CblasConjNoTrans) trans = 0;
    اگر (Trans == CblasConjTrans)   trans = 1;
#وگرنہ_اگر !defined(HEMM)
    اگر (Trans == CblasTrans)       trans = 1;
#ورنہ
    اگر (Trans == CblasConjTrans)   trans = 1;
#اگر_ختم

    info = -1;

    nrowa = دلائل.n;
    اگر (trans & 1) nrowa = دلائل.k;

    اگر (دلائل.ldc < MAX(1,دلائل.n)) info = 10;
    اگر (دلائل.lda < MAX(1,nrowa))  info =  7;
    اگر (دلائل.k < 0)               info =  4;
    اگر (دلائل.n < 0)               info =  3;
    اگر (trans < 0)                info =  2;
    اگر (uplo  < 0)                info =  1;
  }

  اگر (ترتیب == Cblas_صف_Major) {
    اگر (Uplo == Cblas_بالائی) uplo  = 1;
    اگر (Uplo == Cblas_زیریں) uplo  = 0;

    اگر (Trans == CblasNoTrans)     trans = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
    اگر (Trans == CblasTrans)       trans = 0;
    اگر (Trans == CblasConjNoTrans) trans = 1;
    اگر (Trans == CblasConjTrans)   trans = 0;
#وگرنہ_اگر !defined(HEMM)
    اگر (Trans == CblasTrans)       trans = 0;
#ورنہ
    اگر (Trans == CblasConjTrans)   trans = 0;
#اگر_ختم

    info = -1;

    nrowa = دلائل.n;
    اگر (trans & 1) nrowa = دلائل.k;

    اگر (دلائل.ldc < MAX(1,دلائل.n)) info = 10;
    اگر (دلائل.lda < MAX(1,nrowa))  info =  7;
    اگر (دلائل.k < 0)               info =  4;
    اگر (دلائل.n < 0)               info =  3;
    اگر (trans < 0)                info =  2;
    اگر (uplo  < 0)                info =  1;
  }

  اگر (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

  اگر (دلائل.n == 0) واپس;

#اگر_تعریف_شدہ DYNAMIC_ARCH
خارجی حرف* gotoblas_corename(خالی);
#اگر_ختم
	
#اگر !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#اگر defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_براہ_راست)||defined(DYNAMIC_ARCH))
#اگر defined(DYNAMIC_ARCH)
اگر (strcmp(gotoblas_corename(), "armv9sme") == 0
#اگر defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#اگر_ختم
)
#اگر_ختم
  اگر (ترتیب == Cblas_صف_Major && n == ldc) {
    اگر (Trans == CblasNoTrans && k == lda) {
      (Uplo == Cblas_بالائی ? SSYRK_براہ_راست_ALPHA_BETA_UN : SSYRK_براہ_راست_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      واپس;
    } ورنہ اگر (Trans == CblasTrans && n == lda){
      (Uplo == Cblas_بالائی ? SSYRK_براہ_راست_ALPHA_BETA_UT : SSYRK_براہ_راست_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      واپس;
    }
  }
#اگر_ختم
#اگر_ختم

#اگر_ختم


  IDEBUG_آغاز;

  FUNCTION_PROFILE_آغاز();

  بفر = (FLOAT *)blas_حافظہ_alloc(0);
  اگر (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    واپس;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_انحراف_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_انحراف_B);

#اگر_تعریف_شدہ SMP
#اگر_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3
  اگر (!trans){
    منوال |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ورنہ {
    منوال |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  منوال |= (uplo  << BLAS_UPLO_سرکاؤ);
#اگر_ختم

  دلائل.common = تہی;

  NNK = (دوہرا)(دلائل.n+1)*(دوہرا)دلائل.n*(دوہرا)دلائل.k;
  اگر (NNK <= (SMP_حد_MIN * GEMM_MULTITHREAD_حد)) {
  دلائل.nthreads = 1;
  } ورنہ {
  دلائل.nthreads = عدد_cpu_avail(3);
  }

  اگر (دلائل.nthreads == 1) {
#اگر_ختم

    (syrk[(uplo << 1) | trans ])(&دلائل, تہی, تہی, sa, sb, 0);

#اگر_تعریف_شدہ SMP
  } ورنہ {

#اگر_غیر_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&دلائل, تہی, تہی, sa, sb, 0);

#ورنہ

    syrk_thread(منوال, &دلائل, تہی, تہی, syrk[(uplo << 1) | trans ], sa, sb, دلائل.nthreads);

#اگر_ختم

  }
#اگر_ختم

 blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_انتہا(COMPSIZE * COMPSIZE, دلائل.n * دلائل.k + دلائل.n * دلائل.n / 2, دلائل.n * دلائل.n * دلائل.k);

  IDEBUG_انتہا;

  واپس;
}
