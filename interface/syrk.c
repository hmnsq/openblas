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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <ctype.h>
#અન્તર્ભાવ "common.h"
#ifdef FUNCTION_PROFILE
#અન્તર્ભાવ "functable.h"
#endif

#ifndef COMPLEX
#define SMP_THRESHOLD_MIN 109944.
#ifdef XDOUBLE
#define ERROR_NAME "QSYRK "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DSYRK "
#નહીંતર
#define ERROR_NAME "SSYRK "
#endif
#નહીંતર
#define SMP_THRESHOLD_MIN 14824.
#ifndef HEMM
#ifdef XDOUBLE
#define ERROR_NAME "XSYRK "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZSYRK "
#નહીંતર
#define ERROR_NAME "CSYRK "
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "XHERK "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZHERK "
#નહીંતર
#define ERROR_NAME "CHERK "
#endif
#endif
#endif

#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

સ્થિર પૂર્ણાંક (*syrk[])(blas_arg_t *, BLASLONG *, BLASLONG *, FLOAT *, FLOAT *, BLASLONG) = {
#ifndef HEMM
  SYRK_UN, SYRK_UC, SYRK_LN, SYRK_LC,
#જો defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  SYRK_THREAD_UN, SYRK_THREAD_UC, SYRK_THREAD_LN, SYRK_THREAD_LC,
#endif
#નહીંતર
  HERK_UN, HERK_UC, HERK_LN, HERK_LC,
#જો defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)
  HERK_THREAD_UN, HERK_THREAD_UC, HERK_THREAD_LN, HERK_THREAD_LC,
#endif
#endif
};

OPENBLAS_EXPORT

#ifndef CBLAS

શૂન્ય NAME(અક્ષર *UPLO, અક્ષર *TRANS,
         blasint *N, blasint *K,
         FLOAT *alpha, FLOAT *a, blasint *ldA,
         FLOAT *beta,  FLOAT *c, blasint *ldC){

  અક્ષર uplo_arg  = *UPLO;
  અક્ષર trans_arg = *TRANS;

  blas_arg_t args;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
  દ્વિગુણ NNK;
#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_REAL;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_REAL;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#નહીંતર
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  blasint info;
  પૂર્ણાંક uplo;
  પૂર્ણાંક trans;
  પૂર્ણાંક nrowa;

  PRINT_DEBUG_NAME;

  args.n = *N;
  args.k = *K;

  args.a = (શૂન્ય *)a;
  args.c = (શૂન્ય *)c;

  args.lda = *ldA;
  args.ldc = *ldC;

  args.alpha = (શૂન્ય *)alpha;
  args.beta  = (શૂન્ય *)beta;

  TOUPPER(uplo_arg);
  TOUPPER(trans_arg);

  uplo  = -1;
  trans = -1;

  જો (uplo_arg  == 'U') uplo  = 0;
  જો (uplo_arg  == 'L') uplo  = 1;


#ifndef COMPLEX
  જો (trans_arg == 'N') trans = 0;
  જો (trans_arg == 'T') trans = 1;
  જો (trans_arg == 'C') trans = 1;
#નહીંતર
#ifdef HEMM
  જો (trans_arg == 'N') trans = 0;
  જો (trans_arg == 'C') trans = 1;
#નહીંતર
  જો (trans_arg == 'N') trans = 0;
  જો (trans_arg == 'T') trans = 1;
#endif

#endif

  nrowa = args.n;
  જો (trans & 1) nrowa = args.k;

  info = 0;

  જો (args.ldc < MAX(1,args.n)) info = 10;
  જો (args.lda < MAX(1,nrowa))  info =  7;
  જો (args.k < 0)               info =  4;
  જો (args.n < 0)               info =  3;
  જો (trans < 0)                info =  2;
  જો (uplo  < 0)                info =  1;

  જો (info != 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_UPLO Uplo, ગણના CBLAS_TRANSPOSE Trans,
	   blasint n, blasint k,
#જો !defined(COMPLEX) || defined(HEMM)
	   FLOAT alpha,
#નહીંતર
	   શૂન્ય *valpha,
#endif
#જો !defined(COMPLEX)
	   FLOAT *a, blasint lda,
#નહીંતર
	   શૂન્ય *va, blasint lda,
#endif
#જો !defined(COMPLEX) || defined(HEMM)
	   FLOAT beta,
#નહીંતર
	   શૂન્ય *vbeta,
#endif
#જો !defined(COMPLEX)
	   FLOAT *c, blasint ldc) {
#નહીંતર
	   શૂન્ય *vc, blasint ldc) {
#endif

#ifdef COMPLEX
#જો !defined(HEMM)
  FLOAT* alpha = (FLOAT*) valpha;
  FLOAT* beta = (FLOAT*) vbeta;
#endif
  FLOAT* a = (FLOAT*) va;
  FLOAT* c = (FLOAT*) vc;
#endif

  blas_arg_t args;
  પૂર્ણાંક uplo, trans;
  blasint info, nrowa;

  FLOAT *buffer;
  FLOAT *sa, *sb;

#ifdef SMP
દ્વિગુણ NNK;

#ifdef USE_SIMPLE_THREADED_LEVEL3
#ifndef COMPLEX
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_REAL;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_REAL;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#નહીંતર
#ifdef XDOUBLE
  પૂર્ણાંક mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#અથવાજો defined(DOUBLE)
  પૂર્ણાંક mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#નહીંતર
  પૂર્ણાંક mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
#endif

  PRINT_DEBUG_CNAME;

  args.n = n;
  args.k = k;

  args.a = (શૂન્ય *)a;
  args.c = (શૂન્ય *)c;

  args.lda = lda;
  args.ldc = ldc;

#જો !defined(COMPLEX) || defined(HEMM)
  args.alpha = (શૂન્ય *)&alpha;
  args.beta  = (શૂન્ય *)&beta;
#નહીંતર
  args.alpha = (શૂન્ય *)alpha;
  args.beta  = (શૂન્ય *)beta;
#endif

  trans = -1;
  uplo  = -1;
  info  =  0;

  જો (order == CblasColMajor) {
    જો (Uplo == CblasUpper) uplo  = 0;
    જો (Uplo == CblasLower) uplo  = 1;

    જો (Trans == CblasNoTrans)     trans = 0;
#ifndef COMPLEX
    જો (Trans == CblasTrans)       trans = 1;
    જો (Trans == CblasConjNoTrans) trans = 0;
    જો (Trans == CblasConjTrans)   trans = 1;
#અથવાજો !defined(HEMM)
    જો (Trans == CblasTrans)       trans = 1;
#નહીંતર
    જો (Trans == CblasConjTrans)   trans = 1;
#endif

    info = -1;

    nrowa = args.n;
    જો (trans & 1) nrowa = args.k;

    જો (args.ldc < MAX(1,args.n)) info = 10;
    જો (args.lda < MAX(1,nrowa))  info =  7;
    જો (args.k < 0)               info =  4;
    જો (args.n < 0)               info =  3;
    જો (trans < 0)                info =  2;
    જો (uplo  < 0)                info =  1;
  }

  જો (order == CblasRowMajor) {
    જો (Uplo == CblasUpper) uplo  = 1;
    જો (Uplo == CblasLower) uplo  = 0;

    જો (Trans == CblasNoTrans)     trans = 1;
#ifndef COMPLEX
    જો (Trans == CblasTrans)       trans = 0;
    જો (Trans == CblasConjNoTrans) trans = 1;
    જો (Trans == CblasConjTrans)   trans = 0;
#અથવાજો !defined(HEMM)
    જો (Trans == CblasTrans)       trans = 0;
#નહીંતર
    જો (Trans == CblasConjTrans)   trans = 0;
#endif

    info = -1;

    nrowa = args.n;
    જો (trans & 1) nrowa = args.k;

    જો (args.ldc < MAX(1,args.n)) info = 10;
    જો (args.lda < MAX(1,nrowa))  info =  7;
    જો (args.k < 0)               info =  4;
    જો (args.n < 0)               info =  3;
    જો (trans < 0)                info =  2;
    જો (uplo  < 0)                info =  1;
  }

  જો (info >= 0) {
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

  જો (args.n == 0) પાછા;

#ifdef DYNAMIC_ARCH
બાહ્ય અક્ષર* gotoblas_corename(શૂન્ય);
#endif
	
#જો !defined(COMPLEX) && !defined(DOUBLE) && !defined(BFLOAT16)  && !defined(HFLOAT16)
#જો defined(ARCH_ARM64) && (defined(USE_SSYRK_KERNEL_DIRECT)||defined(DYNAMIC_ARCH))
#જો defined(DYNAMIC_ARCH)
જો (strcmp(gotoblas_corename(), "armv9sme") == 0
#જો defined(__clang__)
 || strcmp(gotoblas_corename(), "vortexm4") == 0
#endif
)
#endif
  જો (order == CblasRowMajor && n == ldc) {
    જો (Trans == CblasNoTrans && k == lda) {
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UN : SSYRK_DIRECT_ALPHA_BETA_LN)(n, k, alpha, a, lda, beta, c, ldc);
      પાછા;
    } નહીંતર જો (Trans == CblasTrans && n == lda){
      (Uplo == CblasUpper ? SSYRK_DIRECT_ALPHA_BETA_UT : SSYRK_DIRECT_ALPHA_BETA_LT)(n, k, alpha, a, lda, beta, c, ldc);
      પાછા;
    }
  }
#endif
#endif

#endif


  IDEBUG_START;

  FUNCTION_PROFILE_START();

  buffer = (FLOAT *)blas_memory_alloc(0);
  જો (!buffer) {
    info = -999;
    BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
    પાછા;
  }

  sa = (FLOAT *)((BLASLONG)buffer + GEMM_OFFSET_A);
  sb = (FLOAT *)(((BLASLONG)sa + ((GEMM_P * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)) + GEMM_OFFSET_B);

#ifdef SMP
#ifdef USE_SIMPLE_THREADED_LEVEL3
  જો (!trans){
    mode |= (BLAS_TRANSA_N | BLAS_TRANSB_T);
  } નહીંતર {
    mode |= (BLAS_TRANSA_T | BLAS_TRANSB_N);
  }
  mode |= (uplo  << BLAS_UPLO_SHIFT);
#endif

  args.common = ખાલી;

  NNK = (દ્વિગુણ)(args.n+1)*(દ્વિગુણ)args.n*(દ્વિગુણ)args.k;
  જો (NNK <= (SMP_THRESHOLD_MIN * GEMM_MULTITHREAD_THRESHOLD)) {
  args.nthreads = 1;
  } નહીંતર {
  args.nthreads = અંક_cpu_avail(3);
  }

  જો (args.nthreads == 1) {
#endif

    (syrk[(uplo << 1) | trans ])(&args, ખાલી, ખાલી, sa, sb, 0);

#ifdef SMP
  } નહીંતર {

#ifndef USE_SIMPLE_THREADED_LEVEL3

    (syrk[4 | (uplo << 1) | trans ])(&args, ખાલી, ખાલી, sa, sb, 0);

#નહીંતર

    syrk_thread(mode, &args, ખાલી, ખાલી, syrk[(uplo << 1) | trans ], sa, sb, args.nthreads);

#endif

  }
#endif

 blas_memory_free(buffer);

  FUNCTION_PROFILE_END(COMPSIZE * COMPSIZE, args.n * args.k + args.n * args.n / 2, args.n * args.n * args.k);

  IDEBUG_END;

  પાછા;
}
