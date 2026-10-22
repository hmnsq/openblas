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

#அந்தர்பாவ <stdio.h>
#அந்தர்பாவ <ctype.h>
#அந்தர்பாவ "common.h"
#ifdef FUNCTION_PROFILE
#அந்தர்பாவ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 109944.
#ifdef XDOUBLE
#define ERROR_NAME "QSYRK "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "DSYRK "
#இல்லையெனில்
#define ERROR_NAME "SSYRK "
#endif
#இல்லையெனில்
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#இல்லையெனில்
#define ERROR_NAME "CSYRK "
#endif
#இல்லையெனில்
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#அல்லதெனில் defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#இல்லையெனில்
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

நிலையான முழு (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#இல்லையெனில்
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#எனில் defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

வெறுமை NAME(எழுத்து *UPLO, எழுத்து *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  எழுத்து uplo_arg  = *UPLO;
  எழுத்து trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  இரட்டை NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  முழு uplo;
  முழு trans;
  முழு nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (வெறுமை *)a;
  args.c = (வெறுமை *)c;

  args.lda = *ldA;
  args.ldc = *ldC;

  args.alpha = (வெறுமை *)alpha;
  args.beta  = (வெறுமை *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  எனில் (uplo_arg  == 'U') uplo  = 0;
  எனில் (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'T') trans = 1;
  எனில் (trans_arg == 'C') trans = 1;
#இல்லையெனில்
#ifdef HEMM
  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'C') trans = 1;
#இல்லையெனில்
  எனில் (trans_arg == 'N') trans = 0;
  எனில் (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = args.n;
  எனில் (trans & 1) nrowa = args.k;

  info = 0;

  எனில் (args.ldc < MAX(1,args.n)) info = 10;
  எனில் (args.lda < MAX(1,nrowa))  info =  7;
  எனில் (args.k < 0)               info =  4;
  எனில் (args.n < 0)               info =  3;
  எனில் (trans < 0)                info =  2;
  எனில் (uplo  < 0)                info =  1;

  எனில் (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

#இல்லையெனில்

வெறுமை CNAME(எண்ணிக்கை CBLAS_ORDER order, எண்ணிக்கை CBLAS_UPLO Uplo, எண்ணிக்கை CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#எனில் !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#இல்லையெனில்
	   வெறுமை *valpha,
#endif
#எனில் !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#இல்லையெனில்
	   வெறுமை *va, blasint lda,
#endif
#எனில் !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#இல்லையெனில்
	   வெறுமை *vbeta,
#endif
#எனில் !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#இல்லையெனில்
	   வெறுமை *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#எனில் !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  முழு uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
இரட்டை NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_REAL;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_REAL;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#இல்லையெனில்
#ifdef XDOUBLE
  முழு mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#அல்லதெனில் defined(DOUBLE)
  முழு mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#இல்லையெனில்
  முழு mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (வெறுமை *)a;
  args.c = (வெறுமை *)c;

  args.lda = lda;
  args.ldc = ldc;

#எனில் !defined(COMPLEX) || defined(HEMM)
  args.alpha = (வெறுமை *)&alpha;
  args.beta  = (வெறுமை *)&beta;
#இல்லையெனில்
  args.alpha = (வெறுமை *)alpha;
  args.beta  = (வெறுமை *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  எனில் (order == CblasColMajor) {
    எனில் (Uplo == CblasUpper) uplo  = 0;
    எனில் (Uplo == CblasLower) uplo  = 1;

    எனில் (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    எனில் (Trans == CblasTrans)       trans = 1;
    எனில் (Trans == CblasConjNoTrans) trans = 0;
    எனில் (Trans == CblasConjTrans)   trans = 1;
#அல்லதெனில் !defined(HEMM)
    எனில் (Trans == CblasTrans)       trans = 1;
#இல்லையெனில்
    எனில் (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    எனில் (trans & 1) nrowa = args.k;

    எனில் (args.ldc < MAX(1,args.n)) info = 10;
    எனில் (args.lda < MAX(1,nrowa))  info =  7;
    எனில் (args.k < 0)               info =  4;
    எனில் (args.n < 0)               info =  3;
    எனில் (trans < 0)                info =  2;
    எனில் (uplo  < 0)                info =  1;
  }

  எனில் (order == CblasRowMajor) {
    எனில் (Uplo == CblasUpper) uplo  = 1;
    எனில் (Uplo == CblasLower) uplo  = 0;

    எனில் (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    எனில் (Trans == CblasTrans)       trans = 0;
    எனில் (Trans == CblasConjNoTrans) trans = 1;
    எனில் (Trans == CblasConjTrans)   trans = 0;
#அல்லதெனில் !defined(HEMM)
    எனில் (Trans == CblasTrans)       trans = 0;
#இல்லையெனில்
    எனில் (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    எனில் (trans & 1) nrowa = args.k;

    எனில் (args.ldc < MAX(1,args.n)) info = 10;
    எனில் (args.lda < MAX(1,nrowa))  info =  7;
    எனில் (args.k < 0)               info =  4;
    எனில் (args.n < 0)               info =  3;
    எனில் (trans < 0)                info =  2;
    எனில் (uplo  < 0)                info =  1;
  }

  எனில் (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

  எனில் (args.n == 0) திரும்பு;

#ifdef DYNAMIC_ARCH
வெளி எழுத்து* gotoblas_corename(வெறுமை);
#endif
	
#எனில் !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#எனில் defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#எனில் defined(DYNAMIC_ARCH)
எனில் (strcmp(gotoblas_corename(), "armv9sme") == 0
#எனில் defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  எனில் (order == CblasRowMajor && n == ldc) {
    எனில் (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      திரும்பு;
    } இல்லையெனில் எனில் (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      திரும்பு;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  எனில் (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, அளவு(ERROR_NAME));
    திரும்பு;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  எனில் (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } இல்லையெனில் {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = வெற்று;

  NNK = (இரட்டை)(args.n+1)*(இரட்டை)args.n*(இரட்டை)args.k;
  எனில் (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } இல்லையெனில் {
  args.nthreads = இலக்கம்_cpu_avail(3);
  }

  எனில் (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, வெற்று, வெற்று, sa, sb, 0);

#ifdef SMP
  } இல்லையெனில் {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, வெற்று, வெற்று, sa, sb, 0);

#இல்லையெனில்

    syrk_thread(mode, &args, வெற்று, வெற்று, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  திரும்பு;
}
