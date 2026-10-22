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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <ctype.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ FUNCTION_PROFILE
#ਅਨ੍ਤਰ੍ਭਾਵ "functable.h"
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਹੱਦ_MIN 109944.
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QSYRK "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DSYRK "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SSYRK "
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMP_ਹੱਦ_MIN 14824.
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ HEMM
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XSYRK "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZSYRK "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CSYRK "
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XHERK "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZHERK "
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CHERK "
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ GEMM_MULTITHREAD_ਹੱਦ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_MULTITHREAD_ਹੱਦ 4
#ਜੇ_ਅੰਤ

ਸਥਿਰ ਪੂਰਨ (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#ਜੇ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#ਜੇ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
};

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS

ਖਾਲੀ NAME(ਅੱਖਰ *UPLO, ਅੱਖਰ *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ਅੱਖਰ uplo_arg  = *UPLO;
  ਅੱਖਰ trans_arg = *TRANS;

  blas_arg_t ਦਲੀਲਾਂ;

  FLOAT *ਬਫਰ;
  FLOAT *sa, *sb;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  ਦੁੱਗਣਾ NNK;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ USE_SIMPLE_THREADED_LEVEL3
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  blasint info;
  ਪੂਰਨ uplo;
  ਪੂਰਨ trans;
  ਪੂਰਨ nrowa;

  PRINT_DEBUG_NAME;

  ਦਲੀਲਾਂ.n = *N;
  ਦਲੀਲਾਂ.k = *K;

  ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)a;
  ਦਲੀਲਾਂ.c = (ਖਾਲੀ *)c;

  ਦਲੀਲਾਂ.lda = *ldA;
  ਦਲੀਲਾਂ.ldc = *ldC;

  ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)alpha;
  ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ਜੇ (uplo_arg  == 'U') uplo  = 0;
  ਜੇ (uplo_arg  == 'L') uplo  = 1;


#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
  ਜੇ (trans_arg == 'N') trans = 0;
  ਜੇ (trans_arg == 'T') trans = 1;
  ਜੇ (trans_arg == 'C') trans = 1;
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ HEMM
  ਜੇ (trans_arg == 'N') trans = 0;
  ਜੇ (trans_arg == 'C') trans = 1;
#ਨਹੀਂਤਾਂ
  ਜੇ (trans_arg == 'N') trans = 0;
  ਜੇ (trans_arg == 'T') trans = 1;
#ਜੇ_ਅੰਤ

#ਜੇ_ਅੰਤ

  nrowa = ਦਲੀਲਾਂ.n;
  ਜੇ (trans & 1) nrowa = ਦਲੀਲਾਂ.k;

  info = 0;

  ਜੇ (ਦਲੀਲਾਂ.ldc < MAX(1,ਦਲੀਲਾਂ.n)) info = 10;
  ਜੇ (ਦਲੀਲਾਂ.lda < MAX(1,nrowa))  info =  7;
  ਜੇ (ਦਲੀਲਾਂ.k < 0)               info =  4;
  ਜੇ (ਦਲੀਲਾਂ.n < 0)               info =  3;
  ਜੇ (trans < 0)                info =  2;
  ਜੇ (uplo  < 0)                info =  1;

  ਜੇ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_UPLO Uplo, ਗਿਣਤੀ CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ਜੇ !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ *valpha,
#ਜੇ_ਅੰਤ
#ਜੇ !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ *va, blasint lda,
#ਜੇ_ਅੰਤ
#ਜੇ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ *vbeta,
#ਜੇ_ਅੰਤ
#ਜੇ !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ *vc, blasint ldc) {
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#ਜੇ_ਅੰਤ
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#ਜੇ_ਅੰਤ

  blas_arg_t ਦਲੀਲਾਂ;
  ਪੂਰਨ uplo, trans;
  blasint info, nrowa;

  FLOAT *ਬਫਰ;
  FLOAT *sa, *sb;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
ਦੁੱਗਣਾ NNK;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ USE_SIMPLE_THREADED_LEVEL3
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
  ਪੂਰਨ ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  ਦਲੀਲਾਂ.n = n;
  ਦਲੀਲਾਂ.k = k;

  ਦਲੀਲਾਂ.a = (ਖਾਲੀ *)a;
  ਦਲੀਲਾਂ.c = (ਖਾਲੀ *)c;

  ਦਲੀਲਾਂ.lda = lda;
  ਦਲੀਲਾਂ.ldc = ldc;

#ਜੇ !defined(COMPLEX) || defined(HEMM)
  ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)&alpha;
  ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)&beta;
#ਨਹੀਂਤਾਂ
  ਦਲੀਲਾਂ.alpha = (ਖਾਲੀ *)alpha;
  ਦਲੀਲਾਂ.beta  = (ਖਾਲੀ *)beta;
#ਜੇ_ਅੰਤ

  trans = -1;
  uplo  = -1;
  info  =  0;

  ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo  = 0;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo  = 1;

    ਜੇ (Trans == CblasNoTrans)     trans = 0;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (Trans == CblasTrans)       trans = 1;
    ਜੇ (Trans == CblasConjNoTrans) trans = 0;
    ਜੇ (Trans == CblasConjTrans)   trans = 1;
#ਜਾਂਜੇ !defined(HEMM)
    ਜੇ (Trans == CblasTrans)       trans = 1;
#ਨਹੀਂਤਾਂ
    ਜੇ (Trans == CblasConjTrans)   trans = 1;
#ਜੇ_ਅੰਤ

    info = -1;

    nrowa = ਦਲੀਲਾਂ.n;
    ਜੇ (trans & 1) nrowa = ਦਲੀਲਾਂ.k;

    ਜੇ (ਦਲੀਲਾਂ.ldc < MAX(1,ਦਲੀਲਾਂ.n)) info = 10;
    ਜੇ (ਦਲੀਲਾਂ.lda < MAX(1,nrowa))  info =  7;
    ਜੇ (ਦਲੀਲਾਂ.k < 0)               info =  4;
    ਜੇ (ਦਲੀਲਾਂ.n < 0)               info =  3;
    ਜੇ (trans < 0)                info =  2;
    ਜੇ (uplo  < 0)                info =  1;
  }

  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
    ਜੇ (Uplo == Cblas_ਉੱਤਲਾ) uplo  = 1;
    ਜੇ (Uplo == Cblas_ਹੇਠਲਾ) uplo  = 0;

    ਜੇ (Trans == CblasNoTrans)     trans = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
    ਜੇ (Trans == CblasTrans)       trans = 0;
    ਜੇ (Trans == CblasConjNoTrans) trans = 1;
    ਜੇ (Trans == CblasConjTrans)   trans = 0;
#ਜਾਂਜੇ !defined(HEMM)
    ਜੇ (Trans == CblasTrans)       trans = 0;
#ਨਹੀਂਤਾਂ
    ਜੇ (Trans == CblasConjTrans)   trans = 0;
#ਜੇ_ਅੰਤ

    info = -1;

    nrowa = ਦਲੀਲਾਂ.n;
    ਜੇ (trans & 1) nrowa = ਦਲੀਲਾਂ.k;

    ਜੇ (ਦਲੀਲਾਂ.ldc < MAX(1,ਦਲੀਲਾਂ.n)) info = 10;
    ਜੇ (ਦਲੀਲਾਂ.lda < MAX(1,nrowa))  info =  7;
    ਜੇ (ਦਲੀਲਾਂ.k < 0)               info =  4;
    ਜੇ (ਦਲੀਲਾਂ.n < 0)               info =  3;
    ਜੇ (trans < 0)                info =  2;
    ਜੇ (uplo  < 0)                info =  1;
  }

  ਜੇ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

  ਜੇ (ਦਲੀਲਾਂ.n == 0) ਵਾਪਸ;

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ DYNAMIC_ARCH
ਬਾਹਰੀ ਅੱਖਰ* gotoblas_corename(ਖਾਲੀ);
#ਜੇ_ਅੰਤ
	
#ਜੇ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ਜੇ defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_ਸਿੱਧਾ)||defined(DYNAMIC_ARCH))
#ਜੇ defined(DYNAMIC_ARCH)
ਜੇ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ਜੇ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#ਜੇ_ਅੰਤ
)
#ਜੇ_ਅੰਤ
  ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major && n == ldc) {
    ਜੇ (Trans == CblasNoTrans && k == lda) {
      (Uplo == Cblas_ਉੱਤਲਾ ? SSYRK_ਸਿੱਧਾ_ALPHA_BETA_UN : SSYRK_ਸਿੱਧਾ_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      ਵਾਪਸ;
    } ਨਹੀਂਤਾਂ ਜੇ (Trans == CblasTrans && n == lda){
      (Uplo == Cblas_ਉੱਤਲਾ ? SSYRK_ਸਿੱਧਾ_ALPHA_BETA_UT : SSYRK_ਸਿੱਧਾ_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      ਵਾਪਸ;
    }
  }
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

#ਜੇ_ਅੰਤ


  IDEBUG_ਸ਼ੁਰੂ;

  FUNCTION_PROFILE_ਸ਼ੁਰੂ();

  ਬਫਰ = (FLOAT *)blas_ਸਮ੍ਰਿਤੀ_alloc(0);
  ਜੇ (!ਬਫਰ) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
    ਵਾਪਸ;
  }

  sa = (FLOAT *)((BLASLONG)ਬਫਰ + GEMM_ਵਿੱਥ_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * ਮਾਪ + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_ਵਿੱਥ_B);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ USE_SIMPLE_THREADED_LEVEL3
  ਜੇ (!trans){
    ਬਹੁਲਕ |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ਨਹੀਂਤਾਂ {
    ਬਹੁਲਕ |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  ਬਹੁਲਕ |= (uplo  << BLAS_UPLO_ਖਿਸਕਾਓ);
#ਜੇ_ਅੰਤ

  ਦਲੀਲਾਂ.common = ਨਿਰਾ;

  NNK = (ਦੁੱਗਣਾ)(ਦਲੀਲਾਂ.n+1)*(ਦੁੱਗਣਾ)ਦਲੀਲਾਂ.n*(ਦੁੱਗਣਾ)ਦਲੀਲਾਂ.k;
  ਜੇ (NNK <= (SMP_ਹੱਦ_MIN * GEMM_MULTITHREAD_ਹੱਦ)) {
  ਦਲੀਲਾਂ.nthreads = 1;
  } ਨਹੀਂਤਾਂ {
  ਦਲੀਲਾਂ.nthreads = ਅੰਕ_cpu_avail(3);
  }

  ਜੇ (ਦਲੀਲਾਂ.nthreads == 1) {
#ਜੇ_ਅੰਤ

    (syrk[(uplo << 1) | trans ])(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
  } ਨਹੀਂਤਾਂ {

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, sa, sb, 0);

#ਨਹੀਂਤਾਂ

    syrk_thread(ਬਹੁਲਕ, &ਦਲੀਲਾਂ, ਨਿਰਾ, ਨਿਰਾ, syrk[(uplo << 1) | trans ], sa, sb, ਦਲੀਲਾਂ.nthreads);

#ਜੇ_ਅੰਤ

  }
#ਜੇ_ਅੰਤ

 blas_ਸਮ੍ਰਿਤੀ_free(ਬਫਰ);

  FUNCTION_PROFILE_ਸਿਰਾ(COMPSIZE * COMPSIZE, ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.k + ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.n / 2, ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.n * ਦਲੀਲਾਂ.k);

  IDEBUG_ਸਿਰਾ;

  ਵਾਪਸ;
}
