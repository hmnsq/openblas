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
#ifdef FUNCTION_PROFILE
#ضمن "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 109944.
#ifdef XDOUBLE
#define ERROR_NAME "QSYRK "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DSYRK "
#ܐܚܪܢܐ
#define ERROR_NAME "SSYRK "
#endif
#ܐܚܪܢܐ
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#ܐܚܪܢܐ
#define ERROR_NAME "CSYRK "
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#ܐܚܪܢܐ
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

ساكن صحيح (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#ܐܢ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#ܐܚܪܢܐ
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#ܐܢ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

فراغ NAME(حرف *UPLO, حرف *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  حرف uplo_arg  = *UPLO;
  حرف trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  مزدوج NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ܐܚܪܢܐ
  صحيح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ܐܚܪܢܐ
  صحيح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  صحيح uplo;
  صحيح trans;
  صحيح nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (فراغ *)a;
  args.c = (فراغ *)c;

  args.lda = *ldA;
  args.ldc = *ldC;

  args.alpha = (فراغ *)alpha;
  args.beta  = (فراغ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ܐܢ (uplo_arg  == 'U') uplo  = 0;
  ܐܢ (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  ܐܢ (trans_arg == 'N') trans = 0;
  ܐܢ (trans_arg == 'T') trans = 1;
  ܐܢ (trans_arg == 'C') trans = 1;
#ܐܚܪܢܐ
#ifdef HEMM
  ܐܢ (trans_arg == 'N') trans = 0;
  ܐܢ (trans_arg == 'C') trans = 1;
#ܐܚܪܢܐ
  ܐܢ (trans_arg == 'N') trans = 0;
  ܐܢ (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = args.n;
  ܐܢ (trans & 1) nrowa = args.k;

  info = 0;

  ܐܢ (args.ldc < MAX(1,args.n)) info = 10;
  ܐܢ (args.lda < MAX(1,nrowa))  info =  7;
  ܐܢ (args.k < 0)               info =  4;
  ܐܢ (args.n < 0)               info =  3;
  ܐܢ (trans < 0)                info =  2;
  ܐܢ (uplo  < 0)                info =  1;

  ܐܢ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

#ܐܚܪܢܐ

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_UPLO Uplo, تعداد CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ܐܢ !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#ܐܚܪܢܐ
	   فراغ *valpha,
#endif
#ܐܢ !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#ܐܚܪܢܐ
	   فراغ *va, blasint lda,
#endif
#ܐܢ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ܐܚܪܢܐ
	   فراغ *vbeta,
#endif
#ܐܢ !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#ܐܚܪܢܐ
	   فراغ *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#ܐܢ !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  صحيح uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
مزدوج NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ܐܚܪܢܐ
  صحيح mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ܐܚܪܢܐ
#ifdef XDOUBLE
  صحيح mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  صحيح mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ܐܚܪܢܐ
  صحيح mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (فراغ *)a;
  args.c = (فراغ *)c;

  args.lda = lda;
  args.ldc = ldc;

#ܐܢ !defined(COMPLEX) || defined(HEMM)
  args.alpha = (فراغ *)&alpha;
  args.beta  = (فراغ *)&beta;
#ܐܚܪܢܐ
  args.alpha = (فراغ *)alpha;
  args.beta  = (فراغ *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  ܐܢ (order == CblasColMajor) {
    ܐܢ (Uplo == CblasUpper) uplo  = 0;
    ܐܢ (Uplo == CblasLower) uplo  = 1;

    ܐܢ (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    ܐܢ (Trans == CblasTrans)       trans = 1;
    ܐܢ (Trans == CblasConjNoTrans) trans = 0;
    ܐܢ (Trans == CblasConjTrans)   trans = 1;
#وإلاإذا !defined(HEMM)
    ܐܢ (Trans == CblasTrans)       trans = 1;
#ܐܚܪܢܐ
    ܐܢ (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    ܐܢ (trans & 1) nrowa = args.k;

    ܐܢ (args.ldc < MAX(1,args.n)) info = 10;
    ܐܢ (args.lda < MAX(1,nrowa))  info =  7;
    ܐܢ (args.k < 0)               info =  4;
    ܐܢ (args.n < 0)               info =  3;
    ܐܢ (trans < 0)                info =  2;
    ܐܢ (uplo  < 0)                info =  1;
  }

  ܐܢ (order == CblasRowMajor) {
    ܐܢ (Uplo == CblasUpper) uplo  = 1;
    ܐܢ (Uplo == CblasLower) uplo  = 0;

    ܐܢ (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    ܐܢ (Trans == CblasTrans)       trans = 0;
    ܐܢ (Trans == CblasConjNoTrans) trans = 1;
    ܐܢ (Trans == CblasConjTrans)   trans = 0;
#وإلاإذا !defined(HEMM)
    ܐܢ (Trans == CblasTrans)       trans = 0;
#ܐܚܪܢܐ
    ܐܢ (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    ܐܢ (trans & 1) nrowa = args.k;

    ܐܢ (args.ldc < MAX(1,args.n)) info = 10;
    ܐܢ (args.lda < MAX(1,nrowa))  info =  7;
    ܐܢ (args.k < 0)               info =  4;
    ܐܢ (args.n < 0)               info =  3;
    ܐܢ (trans < 0)                info =  2;
    ܐܢ (uplo  < 0)                info =  1;
  }

  ܐܢ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

  ܐܢ (args.n == 0) ܦܢܐ;

#ifdef DYNAMIC_ARCH
خارجي حرف* gotoblas_corename(فراغ);
#endif
	
#ܐܢ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ܐܢ defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ܐܢ defined(DYNAMIC_ARCH)
ܐܢ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ܐܢ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  ܐܢ (order == CblasRowMajor && n == ldc) {
    ܐܢ (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      ܦܢܐ;
    } ܐܚܪܢܐ ܐܢ (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      ܦܢܐ;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ܐܢ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
    ܦܢܐ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  ܐܢ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ܐܚܪܢܐ {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = عدم;

  NNK = (مزدوج)(args.n+1)*(مزدوج)args.n*(مزدوج)args.k;
  ܐܢ (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } ܐܚܪܢܐ {
  args.nthreads = num_cpu_avail(3);
  }

  ܐܢ (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, عدم, عدم, sa, sb, 0);

#ifdef SMP
  } ܐܚܪܢܐ {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, عدم, عدم, sa, sb, 0);

#ܐܚܪܢܐ

    syrk_thread(mode, &args, عدم, عدم, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  ܦܢܐ;
}
