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
#جيڪڏهن_وصف_ٿيل FUNCTION_PROFILE
#ضمن "functable.h"
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#وصف_ڪريو SMP_حد_MIN 109944.
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QSYRK "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DSYRK "
#وإلا
#وصف_ڪريو ERROR_NAME "SSYRK "
#جيڪڏهن_پورو
#وإلا
#وصف_ڪريو SMP_حد_MIN 14824.
#جيڪڏهن_وصف_نه_ٿيل HEMM
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XSYRK "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZSYRK "
#وإلا
#وصف_ڪريو ERROR_NAME "CSYRK "
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XHERK "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZHERK "
#وإلا
#وصف_ڪريو ERROR_NAME "CHERK "
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_وصف_نه_ٿيل GEMM_MULTITHREAD_حد
#وصف_ڪريو GEMM_MULTITHREAD_حد 4
#جيڪڏهن_پورو

ساكن صحيح (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#جيڪڏهن_وصف_نه_ٿيل HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#جيڪڏهن_پورو
#وإلا
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#إذا defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#جيڪڏهن_پورو
#جيڪڏهن_پورو
};

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t دليلون;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جيڪڏهن_وصف_ٿيل SMP
  مزدوج NNK;
#جيڪڏهن_وصف_ٿيل USE_SIMPLE_THREADED_LEVEL3
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

  blasint info;
  صحيح uplo;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  دليلون.n = *N;
  دليلون.k = *K;

  دليلون.a = (فراغ *)a;
  دليلون.c = (فراغ *)c;

  دليلون.lda = *ldA;
  دليلون.ldc = *ldC;

  دليلون.alpha = (فراغ *)alpha;
  دليلون.beta  = (فراغ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  إذا (uplo_arg  == 'U') uplo  = 0;
  إذا (uplo_arg  == 'L') uplo  = 1;


#جيڪڏهن_وصف_نه_ٿيل COMPLEX
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
#جيڪڏهن_وصف_ٿيل HEMM
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'C') trans = 1;
#وإلا
  إذا (trans_arg == 'N') trans = 0;
  إذا (trans_arg == 'T') trans = 1;
#جيڪڏهن_پورو

#جيڪڏهن_پورو

  nrowa = دليلون.n;
  إذا (trans & 1) nrowa = دليلون.k;

  info = 0;

  إذا (دليلون.ldc < MAX(1,دليلون.n)) info = 10;
  إذا (دليلون.lda < MAX(1,nrowa))  info =  7;
  إذا (دليلون.k < 0)               info =  4;
  إذا (دليلون.n < 0)               info =  3;
  إذا (trans < 0)                info =  2;
  إذا (uplo  < 0)                info =  1;

  إذا (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب, تعداد CBLAS_UPLO Uplo, تعداد CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#إذا !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#وإلا
	   فراغ *valpha,
#جيڪڏهن_پورو
#إذا !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#وإلا
	   فراغ *va, blasint lda,
#جيڪڏهن_پورو
#إذا !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#وإلا
	   فراغ *vbeta,
#جيڪڏهن_پورو
#إذا !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#وإلا
	   فراغ *vc, blasint ldc) {
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل COMPLEX
#إذا !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#جيڪڏهن_پورو
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#جيڪڏهن_پورو

  blas_arg_t دليلون;
  صحيح uplo, trans;
  blasint info, nrowa;

  FLOAT *بفر;
  FLOAT *sa, *sb;

#جيڪڏهن_وصف_ٿيل SMP
مزدوج NNK;

#جيڪڏهن_وصف_ٿيل USE_SIMPLE_THREADED_LEVEL3
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
  صحيح گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  صحيح گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

  دليلون.n = n;
  دليلون.k = k;

  دليلون.a = (فراغ *)a;
  دليلون.c = (فراغ *)c;

  دليلون.lda = lda;
  دليلون.ldc = ldc;

#إذا !defined(COMPLEX) || defined(HEMM)
  دليلون.alpha = (فراغ *)&alpha;
  دليلون.beta  = (فراغ *)&beta;
#وإلا
  دليلون.alpha = (فراغ *)alpha;
  دليلون.beta  = (فراغ *)beta;
#جيڪڏهن_پورو

  trans = -1;
  uplo  = -1;
  info  =  0;

  إذا (ترتيب == Cblas_کانو_Major) {
    إذا (Uplo == Cblas_مٿانهون) uplo  = 0;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 1;

    إذا (Trans == CblasNoTrans)     trans = 0;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
    إذا (Trans == CblasTrans)       trans = 1;
    إذا (Trans == CblasConjNoTrans) trans = 0;
    إذا (Trans == CblasConjTrans)   trans = 1;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 1;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 1;
#جيڪڏهن_پورو

    info = -1;

    nrowa = دليلون.n;
    إذا (trans & 1) nrowa = دليلون.k;

    إذا (دليلون.ldc < MAX(1,دليلون.n)) info = 10;
    إذا (دليلون.lda < MAX(1,nrowa))  info =  7;
    إذا (دليلون.k < 0)               info =  4;
    إذا (دليلون.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (ترتيب == Cblas_صف_Major) {
    إذا (Uplo == Cblas_مٿانهون) uplo  = 1;
    إذا (Uplo == Cblas_هيٺاهون) uplo  = 0;

    إذا (Trans == CblasNoTrans)     trans = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
    إذا (Trans == CblasTrans)       trans = 0;
    إذا (Trans == CblasConjNoTrans) trans = 1;
    إذا (Trans == CblasConjTrans)   trans = 0;
#وإلاإذا !defined(HEMM)
    إذا (Trans == CblasTrans)       trans = 0;
#وإلا
    إذا (Trans == CblasConjTrans)   trans = 0;
#جيڪڏهن_پورو

    info = -1;

    nrowa = دليلون.n;
    إذا (trans & 1) nrowa = دليلون.k;

    إذا (دليلون.ldc < MAX(1,دليلون.n)) info = 10;
    إذا (دليلون.lda < MAX(1,nrowa))  info =  7;
    إذا (دليلون.k < 0)               info =  4;
    إذا (دليلون.n < 0)               info =  3;
    إذا (trans < 0)                info =  2;
    إذا (uplo  < 0)                info =  1;
  }

  إذا (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  إذا (دليلون.n == 0) أرجع;

#جيڪڏهن_وصف_ٿيل DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#جيڪڏهن_پورو
	
#إذا !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#إذا defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_سڌو)||defined(DYNAMIC_ARCH))
#إذا defined(DYNAMIC_ARCH)
إذا (strcmp(gotoblas_corename(), "armv9sme") == 0
#إذا defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#جيڪڏهن_پورو
)
#جيڪڏهن_پورو
  إذا (ترتيب == Cblas_صف_Major && n == ldc) {
    إذا (Trans == CblasNoTrans && k == lda) {
      (Uplo == Cblas_مٿانهون ? SSYRK_سڌو_ALPHA_BETA_UN : SSYRK_سڌو_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      أرجع;
    } وإلا إذا (Trans == CblasTrans && n == lda){
      (Uplo == Cblas_مٿانهون ? SSYRK_سڌو_ALPHA_BETA_UT : SSYRK_سڌو_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      أرجع;
    }
  }
#جيڪڏهن_پورو
#جيڪڏهن_پورو

#جيڪڏهن_پورو


  IDEBUG_شروعات;

  FUNCTION_PROFILE_شروعات();

  بفر = (FLOAT *)blas_يادگيري_alloc(0);
  إذا (!بفر) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    أرجع;
  }

  sa = (FLOAT *)((BLASLONG)بفر + GEMM_وٿي_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * سائيز + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_وٿي_B);

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_ٿيل USE_SIMPLE_THREADED_LEVEL3
  إذا (!trans){
    گهڻائي |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } وإلا {
    گهڻائي |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  گهڻائي |= (uplo  << BLAS_UPLO_سرڪايو);
#جيڪڏهن_پورو

  دليلون.common = عدم;

  NNK = (مزدوج)(دليلون.n+1)*(مزدوج)دليلون.n*(مزدوج)دليلون.k;
  إذا (NNK <= (SMP_حد_MIN * GEMM_MULTITHREAD_حد)) {
  دليلون.nthreads = 1;
  } وإلا {
  دليلون.nthreads = انگ_cpu_avail(3);
  }

  إذا (دليلون.nthreads == 1) {
#جيڪڏهن_پورو

    (syrk[(uplo << 1) | trans ])(&دليلون, عدم, عدم, sa, sb, 0);

#جيڪڏهن_وصف_ٿيل SMP
  } وإلا {

#جيڪڏهن_وصف_نه_ٿيل USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&دليلون, عدم, عدم, sa, sb, 0);

#وإلا

    syrk_thread(گهڻائي, &دليلون, عدم, عدم, syrk[(uplo << 1) | trans ], sa, sb, دليلون.nthreads);

#جيڪڏهن_پورو

  }
#جيڪڏهن_پورو

 blas_يادگيري_free(بفر);

  FUNCTION_PROFILE_پڇاڙي(COMPSIZE * COMPSIZE, دليلون.n * دليلون.k + دليلون.n * دليلون.n / 2, دليلون.n * دليلون.n * دليلون.k);

  IDEBUG_پڇاڙي;

  أرجع;
}
