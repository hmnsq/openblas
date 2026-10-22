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

#ضمن <stdio.h>
#ضمن <ctype.h>
#ضمن "common.h"
#جے_تعریف_شدہ FUNCTION_PROFILE
#ضمن "functable.h"
#جے_انت

#جے_بنا_تعریف COMPLEX
#تعریف_کرو SMP_حد_MIN 109944.
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QSYRK "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DSYRK "
#وإلا
#تعریف_کرو ERROR_NAME "SSYRK "
#جے_انت
#وإلا
#تعریف_کرو SMP_حد_MIN 14824.
#جے_بنا_تعریف HEMM
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XSYRK "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZSYRK "
#وإلا
#تعریف_کرو ERROR_NAME "CSYRK "
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XHERK "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZHERK "
#وإلا
#تعریف_کرو ERROR_NAME "CHERK "
#جے_انت
#جے_انت
#جے_انت

#جے_بنا_تعریف GEMM_MULTITHREAD_حد
#تعریف_کرو GEMM_MULTITHREAD_حد 4
#جے_انت

ساكن صحيح (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جے_بنا_تعریف HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#جے_انت
#وإلا
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#جے_انت
#جے_انت
};

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t دلیلاں;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جے_تعریف_شدہ SMP
  مزدوج NNK;
#جے_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت
#جے_انت
#جے_انت
#جے_انت

  blasint info;
  صحيح uplo;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  دلیلاں.n = *N;
  دلیلاں.k = *K;

  دلیلاں.a = (فراغ *)a;
  دلیلاں.c = (فراغ *)c;

  دلیلاں.lda = *ldA;
  دلیلاں.ldc = *ldC;

  دلیلاں.alpha = (فراغ *)alpha;
  دلیلاں.beta  = (فراغ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;


#جے_بنا_تعریف COMPLEX
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
#جے_تعریف_شدہ HEMM
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
#جے_انت

#جے_انت

  nrowa = دلیلاں.n;
  إذا (trans & 1) nrowa = دلیلاں.k;

  info = 0;

  إذا (دلیلاں.ldc < MAX(1,دلیلاں.n)) info = 10;
  إذا (دلیلاں.lda < MAX(1,nrowa))  info =  7;
  إذا (دلیلاں.k < 0)               info =  4;
  إذا (دلیلاں.n < 0)               info =  3;
  إذا (trans < 0)                info =  2;
  إذا (uplo  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_UPLO Uplo, تعداد CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#إذا !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#وإلا
	   فراغ *valpha,
#جے_انت
#إذا !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#وإلا
	   فراغ *va, blasint lda,
#جے_انت
#إذا !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#وإلا
	   فراغ *vbeta,
#جے_انت
#إذا !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#وإلا
	   فراغ *vc, blasint ldc) {
#جے_انت

#جے_تعریف_شدہ COMPLEX
#إذا !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#جے_انت
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#جے_انت

  blas_arg_t دلیلاں;
  صحيح uplo, trans;
  blasint info, nrowa;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جے_تعریف_شدہ SMP
مزدوج NNK;

#جے_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_REAL;
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
  صحيح منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت
#جے_انت
#جے_انت
#جے_انت

  PRINT_DEBUG_CNAME;

  دلیلاں.n = n;
  دلیلاں.k = k;

  دلیلاں.a = (فراغ *)a;
  دلیلاں.c = (فراغ *)c;

  دلیلاں.lda = lda;
  دلیلاں.ldc = ldc;

#إذا !defined(COMPLEX) || defined(HEMM)
  دلیلاں.alpha = (فراغ *)&alpha;
  دلیلاں.beta  = (فراغ *)&beta;
#وإلا
  دلیلاں.alpha = (فراغ *)alpha;
  دلیلاں.beta  = (فراغ *)beta;
#جے_انت

  trans = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتیب == Cblas_خانہ_Major) {
    إذا (Uplo == Cblas_اتلا) uplo  = 0;
    إذا (Uplo == Cblas_ہیٹھلا) uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
#جے_بنا_تعریف COMPLEX
    إذا (Trans == CblasTrans)       trans = 1;
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 1;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 1;
#جے_انت

    info = -1;

    nrowa = دلیلاں.n;
    إذا (trans & 1) nrowa = دلیلاں.k;

    إذا (دلیلاں.ldc < MAX(1,دلیلاں.n)) info = 10;
    إذا (دلیلاں.lda < MAX(1,nrowa))  info =  7;
    إذا (دلیلاں.k < 0)               info =  4;
    إذا (دلیلاں.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (ترتیب == Cblas_صف_Major) {
    إذا (Uplo == Cblas_اتلا) uplo  = 1;
    إذا (Uplo == Cblas_ہیٹھلا) uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 1;
#جے_بنا_تعریف COMPLEX
    إذا (Trans == CblasTrans)       trans = 0;
    إذا (Trans == CblasConjNoTrans) trans = 1;
    إذا (Trans == CblasConjTrans)   trans = 0;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 0;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 0;
#جے_انت

    info = -1;

    nrowa = دلیلاں.n;
    إذا (trans & 1) nrowa = دلیلاں.k;

    إذا (دلیلاں.ldc < MAX(1,دلیلاں.n)) info = 10;
    إذا (دلیلاں.lda < MAX(1,nrowa))  info =  7;
    إذا (دلیلاں.k < 0)               info =  4;
    إذا (دلیلاں.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (دلیلاں.n == 0) أرجع;

#جے_تعریف_شدہ DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#جے_انت
	
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_سدھا)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جے_انت
)
#جے_انت
  إذا (ترتیب == Cblas_صف_Major && n == ldc) {
    إذا (Trans == CblasNoTrans && k == lda) {
      (Uplo == Cblas_اتلا ? SSYRK_سدھا_ALPHA_BETA_UN : SSYRK_سدھا_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      أرجع;
    } وإلا إذا (Trans == CblasTrans && n == lda){
      (Uplo == Cblas_اتلا ? SSYRK_سدھا_ALPHA_BETA_UT : SSYRK_سدھا_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      أرجع;
    }
  }
#جے_انت
#جے_انت

#جے_انت


  IDEBUG_شروع;

  FUNCTION_PROFILE_شروع();

  بفر = (FLOAT *)blas_حافظہ_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وتھ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ناپ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وتھ_B);

#جے_تعریف_شدہ SMP
#جے_تعریف_شدہ USE_SIMPLE_THREADED_LEVEL3
  إذا (!trans){
    منوال |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } وإلا {
    منوال |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  منوال |= (uplo  << BLAS_UPLO_کھسکاؤ);
#جے_انت

  دلیلاں.common = عدم;

  NNK = (مزدوج)(دلیلاں.n+1)*(مزدوج)دلیلاں.n*(مزدوج)دلیلاں.k;
  إذا (NNK <= (SMP_حد_MIN * GEMM_MULTITHREAD_حد)) {
  دلیلاں.nthreads = 1;
  } وإلا {
  دلیلاں.nthreads = عدد_cpu_avail(3);
  }

  إذا (دلیلاں.nthreads == 1) {
#جے_انت

    (syrk[(uplo << 1) | trans ])(&دلیلاں, عدم, عدم, sa, sb, 0);

#جے_تعریف_شدہ SMP
  } وإلا {

#جے_بنا_تعریف USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&دلیلاں, عدم, عدم, sa, sb, 0);

#وإلا

    syrk_thread(منوال, &دلیلاں, عدم, عدم, syrk[(uplo << 1) | trans ], sa, sb, دلیلاں.nthreads);

#جے_انت

  }
#جے_انت

 blas_حافظہ_free(بفر);

  FUNCTION_PROFILE_سرا(COMPSIZE * COMPSIZE, دلیلاں.n * دلیلاں.k + دلیلاں.n * دلیلاں.n / 2, دلیلاں.n * دلیلاں.n * دلیلاں.k);

  IDEBUG_سرا;

  أرجع;
}
