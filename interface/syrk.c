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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <ctype.h>
#ಅನ್ತರ್ಭಾವ "common.h"
#ifdef FUNCTION_PROFILE
#ಅನ್ತರ್ಭಾವ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 109944.
#ifdef XDOUBLE
#define ERROR_NAME "QSYRK "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DSYRK "
#ಅನ್ಯಥಾ
#define ERROR_NAME "SSYRK "
#endif
#ಅನ್ಯಥಾ
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CSYRK "
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#ಅನ್ಯಥಾ
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#ಅನ್ಯಥಾ
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#ಯದಿ defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

ಶೂನ್ಯ NAME(ವರ್ಣ *UPLO, ವರ್ಣ *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  ವರ್ಣ uplo_arg  = *UPLO;
  ವರ್ಣ trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  ದ್ವಿಗುಣ NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  ಪೂರ್ಣಾಙ್ಕ uplo;
  ಪೂರ್ಣಾಙ್ಕ trans;
  ಪೂರ್ಣಾಙ್ಕ nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (ಶೂನ್ಯ *)a;
  args.c = (ಶೂನ್ಯ *)c;

  args.lda = *ldA;
  args.ldc = *ldC;

  args.alpha = (ಶೂನ್ಯ *)alpha;
  args.beta  = (ಶೂನ್ಯ *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  ಯದಿ (uplo_arg  == 'U') uplo  = 0;
  ಯದಿ (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
  ಯದಿ (trans_arg == 'C') trans = 1;
#ಅನ್ಯಥಾ
#ifdef HEMM
  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'C') trans = 1;
#ಅನ್ಯಥಾ
  ಯದಿ (trans_arg == 'N') trans = 0;
  ಯದಿ (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = args.n;
  ಯದಿ (trans & 1) nrowa = args.k;

  info = 0;

  ಯದಿ (args.ldc < MAX(1,args.n)) info = 10;
  ಯದಿ (args.lda < MAX(1,nrowa))  info =  7;
  ಯದಿ (args.k < 0)               info =  4;
  ಯದಿ (args.n < 0)               info =  3;
  ಯದಿ (trans < 0)                info =  2;
  ಯದಿ (uplo  < 0)                info =  1;

  ಯದಿ (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_UPLO Uplo, ಗಣನಾ CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#ಯದಿ !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *valpha,
#endif
#ಯದಿ !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *va, blasint lda,
#endif
#ಯದಿ !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *vbeta,
#endif
#ಯದಿ !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#ಯದಿ !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  ಪೂರ್ಣಾಙ್ಕ uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
ದ್ವಿಗುಣ NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
  ಪೂರ್ಣಾಙ್ಕ mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (ಶೂನ್ಯ *)a;
  args.c = (ಶೂನ್ಯ *)c;

  args.lda = lda;
  args.ldc = ldc;

#ಯದಿ !defined(COMPLEX) || defined(HEMM)
  args.alpha = (ಶೂನ್ಯ *)&alpha;
  args.beta  = (ಶೂನ್ಯ *)&beta;
#ಅನ್ಯಥಾ
  args.alpha = (ಶೂನ್ಯ *)alpha;
  args.beta  = (ಶೂನ್ಯ *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  ಯದಿ (order == CblasColMajor) {
    ಯದಿ (Uplo == CblasUpper) uplo  = 0;
    ಯದಿ (Uplo == CblasLower) uplo  = 1;

    ಯದಿ (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    ಯದಿ (Trans == CblasTrans)       trans = 1;
    ಯದಿ (Trans == CblasConjNoTrans) trans = 0;
    ಯದಿ (Trans == CblasConjTrans)   trans = 1;
#ನೋಚೇದ್ಯದಿ !defined(HEMM)
    ಯದಿ (Trans == CblasTrans)       trans = 1;
#ಅನ್ಯಥಾ
    ಯದಿ (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    ಯದಿ (trans & 1) nrowa = args.k;

    ಯದಿ (args.ldc < MAX(1,args.n)) info = 10;
    ಯದಿ (args.lda < MAX(1,nrowa))  info =  7;
    ಯದಿ (args.k < 0)               info =  4;
    ಯದಿ (args.n < 0)               info =  3;
    ಯದಿ (trans < 0)                info =  2;
    ಯದಿ (uplo  < 0)                info =  1;
  }

  ಯದಿ (order == CblasRowMajor) {
    ಯದಿ (Uplo == CblasUpper) uplo  = 1;
    ಯದಿ (Uplo == CblasLower) uplo  = 0;

    ಯದಿ (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    ಯದಿ (Trans == CblasTrans)       trans = 0;
    ಯದಿ (Trans == CblasConjNoTrans) trans = 1;
    ಯದಿ (Trans == CblasConjTrans)   trans = 0;
#ನೋಚೇದ್ಯದಿ !defined(HEMM)
    ಯದಿ (Trans == CblasTrans)       trans = 0;
#ಅನ್ಯಥಾ
    ಯದಿ (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    ಯದಿ (trans & 1) nrowa = args.k;

    ಯದಿ (args.ldc < MAX(1,args.n)) info = 10;
    ಯದಿ (args.lda < MAX(1,nrowa))  info =  7;
    ಯದಿ (args.k < 0)               info =  4;
    ಯದಿ (args.n < 0)               info =  3;
    ಯದಿ (trans < 0)                info =  2;
    ಯದಿ (uplo  < 0)                info =  1;
  }

  ಯದಿ (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

  ಯದಿ (args.n == 0) ಪ್ರತ್ಯಾಗಮ;

#ifdef DYNAMIC_ARCH
ಬಾಹ್ಯ ವರ್ಣ* gotoblas_corename(ಶೂನ್ಯ);
#endif
	
#ಯದಿ !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#ಯದಿ defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#ಯದಿ defined(DYNAMIC_ARCH)
ಯದಿ (strcmp(gotoblas_corename(), "armv9sme") == 0
#ಯದಿ defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  ಯದಿ (order == CblasRowMajor && n == ldc) {
    ಯದಿ (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      ಪ್ರತ್ಯಾಗಮ;
    } ಅನ್ಯಥಾ ಯದಿ (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      ಪ್ರತ್ಯಾಗಮ;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  ಯದಿ (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
    ಪ್ರತ್ಯಾಗಮ;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  ಯದಿ (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } ಅನ್ಯಥಾ {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = ರಿಕ್ತ;

  NNK = (ದ್ವಿಗುಣ)(args.n+1)*(ದ್ವಿಗುಣ)args.n*(ದ್ವಿಗುಣ)args.k;
  ಯದಿ (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } ಅನ್ಯಥಾ {
  args.nthreads = ಅಂಕೆ_cpu_avail(3);
  }

  ಯದಿ (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ifdef SMP
  } ಅನ್ಯಥಾ {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, ರಿಕ್ತ, ರಿಕ್ತ, sa, sb, 0);

#ಅನ್ಯಥಾ

    syrk_thread(mode, &args, ರಿಕ್ತ, ರಿಕ್ತ, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  ಪ್ರತ್ಯಾಗಮ;
}
