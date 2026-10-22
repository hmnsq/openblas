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

#অন্তর্ভাব <stdio.h>
#অন্তর্ভাব <ctype.h>
#অন্তর্ভাব "common.h"
#ifdef FUNCTION_PROFILE
#অন্তর্ভাব "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 109944.
#ifdef XDOUBLE
#define ERROR_NAME "QSYRK "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DSYRK "
#অন্যথা
#define ERROR_NAME "SSYRK "
#endif
#অন্যথা
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#অন্যথা
#define ERROR_NAME "CSYRK "
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#অন্যথা
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

স্থির পূর্ণাঙ্ক (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#অন্যথা
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#যদি defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

শূন্য NAME(বর্ণ *UPLO, বর্ণ *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  বর্ণ uplo_arg  = *UPLO;
  বর্ণ trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  দ্বিগুণ NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  পূর্ণাঙ্ক uplo;
  পূর্ণাঙ্ক trans;
  পূর্ণাঙ্ক nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (শূন্য *)a;
  args.c = (শূন্য *)c;

  args.lda = *ldA;
  args.ldc = *ldC;

  args.alpha = (শূন্য *)alpha;
  args.beta  = (শূন্য *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  যদি (uplo_arg  == 'U') uplo  = 0;
  যদি (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  যদি (trans_arg == 'N') trans = 0;
  যদি (trans_arg == 'T') trans = 1;
  যদি (trans_arg == 'C') trans = 1;
#অন্যথা
#ifdef HEMM
  যদি (trans_arg == 'N') trans = 0;
  যদি (trans_arg == 'C') trans = 1;
#অন্যথা
  যদি (trans_arg == 'N') trans = 0;
  যদি (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = args.n;
  যদি (trans & 1) nrowa = args.k;

  info = 0;

  যদি (args.ldc < MAX(1,args.n)) info = 10;
  যদি (args.lda < MAX(1,nrowa))  info =  7;
  যদি (args.k < 0)               info =  4;
  যদি (args.n < 0)               info =  3;
  যদি (trans < 0)                info =  2;
  যদি (uplo  < 0)                info =  1;

  যদি (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_UPLO Uplo, গণনা CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#যদি !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#অন্যথা
	   শূন্য *valpha,
#endif
#যদি !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#অন্যথা
	   শূন্য *va, blasint lda,
#endif
#যদি !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#অন্যথা
	   শূন্য *vbeta,
#endif
#যদি !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#অন্যথা
	   শূন্য *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#যদি !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  পূর্ণাঙ্ক uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
দ্বিগুণ NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
  পূর্ণাঙ্ক mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
  পূর্ণাঙ্ক mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
  পূর্ণাঙ্ক mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (শূন্য *)a;
  args.c = (শূন্য *)c;

  args.lda = lda;
  args.ldc = ldc;

#যদি !defined(COMPLEX) || defined(HEMM)
  args.alpha = (শূন্য *)&alpha;
  args.beta  = (শূন্য *)&beta;
#অন্যথা
  args.alpha = (শূন্য *)alpha;
  args.beta  = (শূন্য *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  যদি (order == CblasColMajor) {
    যদি (Uplo == CblasUpper) uplo  = 0;
    যদি (Uplo == CblasLower) uplo  = 1;

    যদি (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    যদি (Trans == CblasTrans)       trans = 1;
    যদি (Trans == CblasConjNoTrans) trans = 0;
    যদি (Trans == CblasConjTrans)   trans = 1;
#নোচেদ্যদি !defined(HEMM)
    যদি (Trans == CblasTrans)       trans = 1;
#অন্যথা
    যদি (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    যদি (trans & 1) nrowa = args.k;

    যদি (args.ldc < MAX(1,args.n)) info = 10;
    যদি (args.lda < MAX(1,nrowa))  info =  7;
    যদি (args.k < 0)               info =  4;
    যদি (args.n < 0)               info =  3;
    যদি (trans < 0)                info =  2;
    যদি (uplo  < 0)                info =  1;
  }

  যদি (order == CblasRowMajor) {
    যদি (Uplo == CblasUpper) uplo  = 1;
    যদি (Uplo == CblasLower) uplo  = 0;

    যদি (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    যদি (Trans == CblasTrans)       trans = 0;
    যদি (Trans == CblasConjNoTrans) trans = 1;
    যদি (Trans == CblasConjTrans)   trans = 0;
#নোচেদ্যদি !defined(HEMM)
    যদি (Trans == CblasTrans)       trans = 0;
#অন্যথা
    যদি (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    যদি (trans & 1) nrowa = args.k;

    যদি (args.ldc < MAX(1,args.n)) info = 10;
    যদি (args.lda < MAX(1,nrowa))  info =  7;
    যদি (args.k < 0)               info =  4;
    যদি (args.n < 0)               info =  3;
    যদি (trans < 0)                info =  2;
    যদি (uplo  < 0)                info =  1;
  }

  যদি (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

  যদি (args.n == 0) প্রত্যাগম;

#ifdef DYNAMIC_ARCH
বাহ্য বর্ণ* gotoblas_corename(শূন্য);
#endif
	
#যদি !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#যদি defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#যদি defined(DYNAMIC_ARCH)
যদি (strcmp(gotoblas_corename(), "armv9sme") == 0
#যদি defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  যদি (order == CblasRowMajor && n == ldc) {
    যদি (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      প্রত্যাগম;
    } অন্যথা যদি (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      প্রত্যাগম;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  যদি (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
    প্রত্যাগম;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  যদি (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } অন্যথা {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = রিক্ত;

  NNK = (দ্বিগুণ)(args.n+1)*(দ্বিগুণ)args.n*(দ্বিগুণ)args.k;
  যদি (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } অন্যথা {
  args.nthreads = অঙ্ক_cpu_avail(3);
  }

  যদি (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, রিক্ত, রিক্ত, sa, sb, 0);

#ifdef SMP
  } অন্যথা {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, রিক্ত, রিক্ত, sa, sb, 0);

#অন্যথা

    syrk_thread(mode, &args, রিক্ত, রিক্ত, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  প্রত্যাগম;
}
